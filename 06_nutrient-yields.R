## This script analyzes the detailed catch data with the McClanahan et al. biomass model



##### Load Data #####

## Load trip data
trip <- read.csv("data/clean-data/03_TripData.csv")

## Load site data
site <- read.csv("data/raw-data/landing-sites.csv")

# Make fishing restrictions integers
site$Management <- site$Management %>%
  str_replace_all("Fished", "1") %>%
  str_replace_all("Restricted", "2") %>%
  str_replace_all("Unfished Low compliance", "3") %>%
  str_replace_all("Unfished High compliance", "4") %>%
  as.integer()

# Make stock status integers
site$overfished_msy_underfished <- site$overfished_msy_underfished %>%
  str_replace_all("Overfished", "1") %>%
  str_replace_all("Sustainable", "2") %>%
  str_replace_all("Underfished", "3") %>%
  as.integer()

# Rename columns
site <- site %>%
  select(-c(country, source, ID, Country, N_cells, Ecoregion, geo_source, lost_yield)) %>%
  rename(management = Management,
    current.velocity_m.s = Current.velocity,
    salinity_ppt = Salinity,
    travel.popcenter_hrs = Travel.time.to.nearest.population..hrs.,
    travel.market_hrs = Travel.time.to.nearest.market..hrs.,
    reef.area_km2 = area_sqkm,
    yield_tons.km2 = Yiedl_tons.km2,
    B_Bmsy = B_BMSY,
    stock.status = overfished_msy_underfished,
    SST.median_degC = SST_Median,
    SST.skewness = SST_Skewness,
    SST.kurtosis = SST_Kurtosis,
    SST.rise.rate_degC.yr = SSTRateofRise,
    SST.sd_degC = SST_1,
    SST.bimodality = bimodality,
    global.stress = Globa,
    climate.stress = Clima,
    calcite_mol.m3 = Calcite.concentration,
    ph = PH, 
    DO_mg.l = Dissolved.oxygen,
    diffuse.attenuation.coefficient_m = Diffuse.attenuation.coefficient,
    NPP_mg.m2.day = Net.primary.productivity,
    wave.energy_kW.m = Wave.energy)


## Combine data for analysis

# Aggregate site to one observation per site
temp <- aggregate(. ~ site, data = site, mean)

# Combine data
trip <- left_join(trip, temp, by = "site")

# Make NA in gear to be other
trip$gear[is.na(trip$gear)] <- "other"

# Make character columns factors
trip$country <- as.factor(trip$country)
trip$site <- as.factor(trip$site)
trip$gear <- as.factor(trip$gear)
trip$habitat <- as.factor(trip$habitat)

# Remove rows with missing nutrient values
trip <- trip[!(is.na(trip$nutrients_pdv)),]

# Remove specific nutrient values in favor of cumulative value
trip <- trip %>% select(-(calcium_mg:zinc_pdv.pue))

# Remove additional nutrient values
trip <- trip %>% select(-nutrient_evenness) %>%
  select(-nutrients_pdv.pue)




##### Variable Selection #####


## Add random numbers and remove anything performing worse

# Add random numbers
trip$random <- runif(n = nrow(trip), min = 1, max = 100)

# Data for analysis
trip2 <- select(trip, -c("trip", "date", "latitude", "longitude", "source"))

# Bootstrap sample for model training
trip2.train <- slice_sample(trip2, prop = 1, replace = TRUE)

# Run random forest model
rfm_trip <- randomForest(nutrients_pdv ~ ., data = trip2.train,
  ntree = 200,
  importance = TRUE, proximity = TRUE, na.action = na.roughfix)

# Summary of model
rfm_trip

# 84.99 % var explained on 11/8/23 11:45 (plotted as #01)
# 84.26 % var explained on 11/8/23 11:50 (plotted as #02)
# 82.88 % var explained on 11/8/23 11:52 (plotted as #03)
# 84.2 % var explained on 11/8/23 11:54 (plotted as #04)
# 83.05 % var explained on 11/8/23 11:55 (plotted as #05)

# Get variable importance
ImpData <- as.data.frame(importance(rfm_trip))

# Visualize variable importance
ImpData$Var.Names <- row.names(ImpData)
ImpData$Var.Names <- c("Country", "Site", "Gear", "Landings", "CPUE", "Maturity",
  "Length : Optimum Length", "Length", "K", "Trophic Level",
  "Habitat", "Species Richness", "Management", "Current",
  "Salinity", "Travel to Population", "Travel to Market", "Reef Area", "Yield per Area",
  "B : B[msy]", "Stock Status", "SST Median", "SST Skewness", "SST Kurtosis", "SST Rate of Rise",
  "SST Standard Deviation", "SST Bimodality", "Degree Heating Weeks", "Global Stress",
  "Climate Stress", "PAR max", "Calcite", "ph", "Dissolved Oxygen", "Diffuse Attenuation Coefficient",
  "Net Primary Productivity", "Wave Energy", "Random")
colnames(ImpData)[colnames(ImpData) == "IncNodePurity"] <- "Increase in Node Purity"
ImpData$Var.Names <- fct_reorder(ImpData$Var.Names, ImpData$`%IncMSE`)
ggplot(ImpData, aes(x=Var.Names, y=`%IncMSE`)) +
  geom_segment( aes(x=Var.Names, xend=Var.Names, y=0, yend=`%IncMSE`), color="skyblue") +
  geom_point(aes(size = `Increase in Node Purity`), color="blue", alpha=0.6) +
  theme_light() +
  coord_flip() +
  xlab("") +
  ylab("% Increase in Mean Squared Error") +
  theme_pubr() +
  theme(plot.margin = unit(c(5, 5, 0, 0), "mm"),
    legend.position = "bottom",
    text = element_text(size = 7),
    axis.text.y = element_text(size = 5),
    axis.text.x = element_text(size = 5))

# Save variable importance
ggsave("figures/06_RandomForestAllVariables05.pdf",
  device = "pdf",
  width = 180,
  units = "mm",
  dpi = 300)

# Retain only those variables that performed better than random numbers
trip <- select(trip, c(trip, date, latitude, longitude, source, nutrients_pdv,
  troph, k, habitat, cpue_kg.effort.day, LLmat, length_cm))


## Spearman correlation matrix for numeric variables

# Subset of numeric data
x <- select(trip, troph:length_cm)
x <- select(x, -habitat)

# Correlation matrix
cor.matrix <- cor(x, use = "pairwise.complete.obs", method = "spearman")

# Plot correlations
corrplot(cor.matrix)

# Make correlation matrix a data frame
cor.matrix <- data.frame(cor.matrix)

# Turn all 1's into NA
for(i in 1:nrow(cor.matrix)){
  
  cor.matrix[i, i] <- NA
  
}

# Make all values absolute value
cor.matrix <- abs(cor.matrix)

# Remove one member of all pairs with correlation coefficient above 0.75

# There is nothing to remove. The highest coefficient is
# 0.61 b/w k and troph. Length and L/Lmat are 0.58. Everything
# else is below 0.5.




##### RFM 01 #####

# Data for analysis
trip2 <- select(trip, -c("trip", "date", "latitude", "longitude", "source"))

# For reproducibility
set.seed(1234)

# Impute NAs using random forest proximity
trip2 <- rfImpute(nutrients_pdv ~ ., data = trip2)

# For reproducibility
set.seed(1234)

# Bootstrap sample for model training
trip2.train <- slice_sample(trip2, prop = 1, replace = TRUE)

# The out of bag sample for testing
trip2.test <- setdiff(trip2, trip2.train)

# For reproducibility
set.seed(1234)

# Random forest model with training and testing data
x <- select(trip2.train, -nutrients_pdv)
y <- trip2.train$nutrients_pdv
xtest <- select(trip2.test, -nutrients_pdv)
ytest <- trip2.test$nutrients_pdv
rfm_trip <- randomForest(x = x, y = y, xtest = xtest, ytest = ytest,
  mtry = 4, ntree = 501,
  importance = TRUE, proximity = TRUE, keep.forest = TRUE)

# Summary of model
rfm_trip

# 83.28 train / 68 test % var explained at mtry = 2
# 83.72 train / 68.87 test % var explained at mtry = 3
# 83.89 train / 69.04 test % var explained at mtry = 4
# 83.74 train / 69.08 test % var explained at mtry = 5
# 83.45 train / 68.57 test % var explained at mtry = 6

# We select mtry = 4 because it has the highest performance overall

# Visualize variable importance
ImpData <- as.data.frame(importance(rfm_trip))
ImpData$Var.Names <- row.names(ImpData)
ImpData$Var.Names <- c("Trophic Level", "K", "Habitat", "CPUE", "Maturity", "Length")
ImpData$Var.Names <- fct_reorder(ImpData$Var.Names, ImpData$`%IncMSE`)
colnames(ImpData)[colnames(ImpData) == "IncNodePurity"] <- "Increase in Node Purity"
ggplot(ImpData, aes(x=Var.Names, y=`%IncMSE`)) +
  geom_segment( aes(x=Var.Names, xend=Var.Names, y=0, yend=`%IncMSE`), color="skyblue") +
  geom_point(aes(size = `Increase in Node Purity`), color="blue", alpha=0.6) +
  theme_light() +
  coord_flip() +
  xlab("") +
  ylab("% Increase in Mean Squared Error") +
  theme_pubr() +
  theme(legend.position = "bottom",
    text = element_text(size = 6.5),
    axis.title.x = element_text(vjust = -3))

# # Save plot
# ggsave("figures/06_RFM_VariableImportance.pdf",
#   device = "pdf",
#   width = 114,
#   units = "mm",
#   dpi = 300)




##### Partial Plots #####

## K
pd_k <- bind_rows(partialPlot(rfm_trip, pred.data = x, x.var = "k",
  plot = FALSE, n.pt = 200))
ggplot(pd_k, aes(x = x, y = 100*y)) +
  geom_line() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  coord_cartesian(xlim = c(0, 2), ylim = c(160, 230)) +
  xlab("K") +
  ylab("Nutrient Density (% DRI)") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(plot.margin = unit(c(0.5, 0.5, 1.6, 0.5), "cm"),
    text = element_text(size = 6.5))


## Habitat
pd_habitat <- bind_rows(partialPlot(rfm_trip, pred.data = x, x.var = "habitat",
  plot = FALSE))
pd_habitat$x <- fct_reorder(pd_habitat$x, pd_habitat$y, .desc = TRUE)
ggplot(pd_habitat, aes(x = x, y = 100*y)) +
  geom_bar(stat = "identity") +
  coord_cartesian(ylim = c(160, 230)) +
  xlab("Habitat") +
  ylab("Nutrient Density (% DRI)") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    text = element_text(size = 6.5),
    axis.text.x = element_text(angle = 45, hjust = 1))


## Trophic Level
pd_troph <- bind_rows(partialPlot(rfm_trip, pred.data = x, x.var = "troph",
  plot = FALSE, n.pt = 200))
ggplot(pd_troph, aes(x = x, y = 100*y)) +
  geom_line() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  coord_cartesian(ylim = c(160, 230)) +
  xlab("Trophic Level") +
  ylab("Nutrient Density (% DRI)") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.9, 0.5), "cm"),
    text = element_text(size = 6.5),
    axis.title.x = element_text(margin = unit(c(0.8, 0, 0, 0), "cm")))

## CPUE

# Subset for only CPUE < 25
x2 <- filter(x, cpue_kg.effort.day < 25)

# Plot relationship
pd_cpue <- bind_rows(partialPlot(rfm_trip, pred.data = x2, x.var = "cpue_kg.effort.day",
  plot = FALSE, n.pt = 200))
ggplot(pd_cpue, aes(x = x, y = 100*y)) +
  geom_line() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  coord_cartesian(ylim = c(160, 230)) +
  xlab("Catch per Unit Effort") +
  ylab("Nutrient Density (% DRI)") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    text = element_text(size = 6.5))

# Save image
ggsave("figures/06_RandomForest_CPUE.pdf",
  device = "pdf",
  width = 90,
  units = "mm",
  dpi = 300)


## Maturity

# Subset for only L/Lmat < 10
x2 <- filter(x, LLmat < 10)

# Plot results
pd_LLmat <- bind_rows(partialPlot(rfm_trip, pred.data = x2, x.var = "LLmat",
  plot = FALSE, n.pt = 200))
ggplot(pd_LLmat, aes(x = x, y = 100*y)) +
  geom_line() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  coord_cartesian(ylim = c(160, 230)) +
  #xlab(expression(paste("Maturity ", bgroup("(", frac(L, L[mat]), ")")))) +
  xlab("Maturity") +
  ylab("Nutrient Density (% DRI)") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    text = element_text(size = 6.5))


## Length

# Subset for only length < 100 cm
x2 <- filter(x, length_cm < 100)

# Plot results
pd_length <- bind_rows(partialPlot(rfm_trip, pred.data = x2, x.var = "length_cm",
  plot = FALSE, n.pt = 200))
ggplot(pd_length, aes(x = x, y = 100*y)) +
  geom_line() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  coord_cartesian(ylim = c(160, 230)) +
  xlab("Length (cm)") +
  ylab("Nutrient Density (% DRI)") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    text = element_text(size = 6.5))


## Evaluate partial dependence plots for shape-based significance (Greenwell et al. 2018).

# K = 0.1074631
sd(pd_k$y)

# Habitat = 0.131694
sd(pd_habitat$y)

# Trophic Level = 0.08777557
sd(pd_troph$y)

# CPUE = 0.007590773 - NOT SIGNIFICANT
sd(pd_cpue$y)

# L/Lmat = 0.02026397
sd(pd_LLmat$y)

# Length = 0.1756916
sd(pd_length$y)




##### RFM 02 - Pruned #####

## A new RFM with CPUE removed

# Remove species richness and CPUE
trip2 <- select(trip2, -c("cpue_kg.effort.day"))

# For reproducibility
set.seed(1234)

# Bootstrap sample for model training
trip2.train <- slice_sample(trip2, prop = 1, replace = TRUE)

# The out of bag sample for testing
trip2.test <- setdiff(trip2, trip2.train)

# For reproducibility
set.seed(1234)

# Random forest model with training and testing data
x <- select(trip2.train, -nutrients_pdv)
y <- trip2.train$nutrients_pdv
xtest <- select(trip2.test, -nutrients_pdv)
ytest <- trip2.test$nutrients_pdv
rfm_trip <- randomForest(x = x, y = y, xtest = xtest, ytest = ytest,
  mtry = 3, ntree = 501,
  importance = TRUE, proximity = TRUE, keep.forest = TRUE)

# Summary of model
rfm_trip

# 84.38 train / 66.15 test % var explained at mtry = 2
# 84.72 train / 66.95 test % var explained at mtry = 3
# 84.51 train / 66.92 test % var explained at mtry = 4
# 84.21 train / 66.64 test % var explained at mtry = 5

# We chose mtry = 3 because it has the best fit

# Visualize variable importance
ImpData <- as.data.frame(importance(rfm_trip))
ImpData$Var.Names <- row.names(ImpData)
ImpData$Var.Names <- c("K", "Habitat", "Trophic Level", "Maturity", "Length")
ImpData$Var.Names <- fct_reorder(ImpData$Var.Names, ImpData$`%IncMSE`)
colnames(ImpData)[colnames(ImpData) == "IncNodePurity"] <- "Increase in Node Purity"
a <- ggplot(ImpData, aes(x=Var.Names, y=`%IncMSE`)) +
  geom_segment( aes(x=Var.Names, xend=Var.Names, y=0, yend=`%IncMSE`), color="skyblue") +
  geom_point(aes(size = `Increase in Node Purity`), color="blue", alpha=0.6) +
  theme_light() +
  coord_flip() +
  xlab("") +
  ylab("% Increase in Mean Squared Error") +
  theme_pubr() +
  theme(legend.position = "bottom",
    text = element_text(size = 7),
    axis.text.y = element_text(size = 5),
    axis.text.x = element_text(size = 5))




##### Partial Plots #####

## K
pd_k <- bind_rows(partialPlot(rfm_trip, pred.data = x, x.var = "k",
  plot = FALSE, n.pt = 200))
d <- ggplot(pd_k, aes(x = x, y = 100*y)) +
  geom_line() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  coord_cartesian(xlim = c(0, 2), ylim = c(160, 230)) +
  xlab("K") +
  ylab("Nutrient Density (% DRI)") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    text = element_text(size = 6.5))


## Habitat
pd_habitat <- bind_rows(partialPlot(rfm_trip, pred.data = x, x.var = "habitat",
  plot = FALSE))
pd_habitat$x <- fct_reorder(pd_habitat$x, pd_habitat$y, .desc = TRUE)
b <- ggplot(pd_habitat, aes(x = x, y = 100*y)) +
  geom_bar(stat = "identity") +
  coord_cartesian(ylim = c(160, 230)) +
  xlab("Habitat") +
  ylab("Nutrient Density (% DRI)") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    text = element_text(size = 6.5),
    axis.text.x = element_text(angle = 45, hjust = 1))


## Trophic Level
pd_troph <- bind_rows(partialPlot(rfm_trip, pred.data = x, x.var = "troph",
  plot = FALSE, n.pt = 200))
c <- ggplot(pd_troph, aes(x = x, y = 100*y)) +
  geom_line() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  coord_cartesian(ylim = c(160, 230)) +
  xlab("Trophic Level") +
  ylab("Nutrient Density (% DRI)") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    text = element_text(size = 6.5))


## Maturity

# Subset for only L/Lmat < 10
x2 <- filter(x, LLmat < 10)

# Plot results
pd_LLmat <- bind_rows(partialPlot(rfm_trip, pred.data = x2, x.var = "LLmat",
  plot = FALSE, n.pt = 200))
e <- ggplot(pd_LLmat, aes(x = x, y = 100*y)) +
  geom_line() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  coord_cartesian(ylim = c(160, 230)) +
  #xlab(expression(paste("Maturity ", bgroup("(", frac(L, L[mat]), ")")))) +
  xlab("Maturity") +
  ylab("Nutrient Density (% DRI)") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    text = element_text(size = 6.5))


## Length

# Subset for only length < 100 cm
x2 <- filter(x, length_cm < 100)

# Plot results
pd_length <- bind_rows(partialPlot(rfm_trip, pred.data = x2, x.var = "length_cm",
  plot = FALSE, n.pt = 200))
f <- ggplot(pd_length, aes(x = x, y = 100*y)) +
  geom_line() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  coord_cartesian(ylim = c(160, 230)) +
  xlab("Length (cm)") +
  ylab("Nutrient Density (% DRI)") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    text = element_text(size = 6.5))


## Evaluate partial dependence plots for shape-based significance (Greenwell et al. 2018).

# K = 0.1050185
sd(pd_k$y)

# Habitat = 0.1276161
sd(pd_habitat$y)

# Trophic Level = 0.0865799
sd(pd_troph$y)

# L/Lmat = 0.02130037
sd(pd_LLmat$y)

# Length = 0.1836424
sd(pd_length$y)

## Combine plots for figure
ggarrange(a, b, c, d, e, f,
  ncol = 2, nrow = 3,
  labels = "AUTO",
  font.label = list(size = 7))

# Save image
ggsave("figures/06_RandomForestResults.pdf",
  device = "pdf",
  width = 172, height = 200,
  units = "mm",
  dpi = 300)




##### Re-Load data #####

## Load nutrition requirements
dri <- read.csv("data/clean-data/02_DietaryReferenceIntakes.csv")

## Load trip data
trip <- read.csv("data/clean-data/03_TripData.csv")

## Load site data
site <- read.csv("data/raw-data/landing-sites.csv")

# Make fishing restrictions integers
site$Management <- site$Management %>%
  str_replace_all("Fished", "1") %>%
  str_replace_all("Restricted", "2") %>%
  str_replace_all("Unfished Low compliance", "3") %>%
  str_replace_all("Unfished High compliance", "4") %>%
  as.integer()

# Make stock status integers
site$overfished_msy_underfished <- site$overfished_msy_underfished %>%
  str_replace_all("Overfished", "1") %>%
  str_replace_all("Sustainable", "2") %>%
  str_replace_all("Underfished", "3") %>%
  as.integer()

# Rename columns
site <- site %>%
  select(-c(country, source, ID, Country, N_cells, Ecoregion, geo_source, lost_yield)) %>%
  rename(management = Management,
    current.velocity_m.s = Current.velocity,
    salinity_ppt = Salinity,
    travel.popcenter_hrs = Travel.time.to.nearest.population..hrs.,
    travel.market_hrs = Travel.time.to.nearest.market..hrs.,
    reef.area_km2 = area_sqkm,
    yield_tons.km2 = Yiedl_tons.km2,
    B_Bmsy = B_BMSY,
    stock.status = overfished_msy_underfished,
    SST.median_degC = SST_Median,
    SST.skewness = SST_Skewness,
    SST.kurtosis = SST_Kurtosis,
    SST.rise.rate_degC.yr = SSTRateofRise,
    SST.sd_degC = SST_1,
    SST.bimodality = bimodality,
    global.stress = Globa,
    climate.stress = Clima,
    calcite_mol.m3 = Calcite.concentration,
    ph = PH, 
    DO_mg.l = Dissolved.oxygen,
    diffuse.attenuation.coefficient_m = Diffuse.attenuation.coefficient,
    NPP_mg.m2.day = Net.primary.productivity,
    wave.energy_kW.m = Wave.energy)


## Combine data for analysis

# Aggregate site to one observation per site
temp <- aggregate(. ~ site, data = site, mean)

# Combine data
trip <- left_join(trip, temp, by = "site")

# Make NA in gear to be other
trip$gear[is.na(trip$gear)] <- "other"

# Make character columns factors
trip$country <- as.factor(trip$country)
trip$site <- as.factor(trip$site)
trip$gear <- as.factor(trip$gear)
trip$habitat <- as.factor(trip$habitat)

# Remove rows with missing nutrient values
trip <- trip[!(is.na(trip$nutrients_pdv)),]




##### RFM 03 - Calcium #####

# Select desired variables only
trip2 <- select(trip, c("troph", "k", "habitat", "LLmat", "length_cm", "calcium_mg.per.100g"))

# For reproducibility
set.seed(1234)

# Bootstrap sample for model training
trip2.train <- slice_sample(trip2, prop = 1, replace = TRUE)

# The out of bag sample for testing
trip2.test <- setdiff(trip2, trip2.train)

# For reproducibility
set.seed(1234)

# Random forest model with training and testing data
x <- select(trip2.train, -calcium_mg.per.100g)
y <- trip2.train$calcium_mg.per.100g
xtest <- select(trip2.test, -calcium_mg.per.100g)
ytest <- trip2.test$calcium_mg.per.100g
rfm_trip <- randomForest(x = x, y = y, xtest = xtest, ytest = ytest,
  mtry = 2, ntree = 501,
  importance = TRUE, proximity = TRUE, keep.forest = TRUE)

# Summary of model
rfm_trip

# 84.15 train / 72.02 test % var explained at mtry = 2
# 84.79 train / 70.99 test % var explained at mtry = 3
# 85.13 train / 70.23 test % var explained at mtry = 4
# 85.2 train / 68.54 test % var explained at mtry = 5

# We chose mtry = 2 because it has the best fit for testing data

# Visualize variable importance
ImpData <- as.data.frame(importance(rfm_trip))
ImpData$Var.Names <- row.names(ImpData)
ImpData$Var.Names <- c("Trophic Level", "K", "Habitat",  "Maturity", "Length")
ImpData$Var.Names <- fct_reorder(ImpData$Var.Names, ImpData$`%IncMSE`)
colnames(ImpData)[colnames(ImpData) == "IncNodePurity"] <- "Increase in Node Purity"
a <- ggplot(ImpData, aes(x=Var.Names, y=`%IncMSE`)) +
  geom_segment( aes(x=Var.Names, xend=Var.Names, y=0, yend=`%IncMSE`), color="skyblue") +
  geom_point(aes(size = `Increase in Node Purity`), color="blue", alpha=0.6) +
  theme_light() +
  coord_flip() +
  xlab("") +
  ylab("% Increase in Mean Squared Error") +
  theme_pubr() +
  theme(legend.position = "bottom",
    plot.margin = unit(c(0.5, 1, 0.5, 0), "cm"),
    text = element_text(size = 7),
    axis.text.y = element_text(size = 5),
    axis.text.x = element_text(size = 5))




##### Partial Plots #####

## K
pd_k <- bind_rows(partialPlot(rfm_trip, pred.data = x, x.var = "k",
  plot = FALSE, n.pt = 200))
c <- ggplot(pd_k, aes(x = x, y = y)) +
  geom_line() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  # geom_hline(yintercept = dri$dri[dri$nutrients == "calcium_mg"],
  #   linetype = "dashed") +
  coord_cartesian(xlim = c(0, 2), ylim = c(30, 140)) +
  xlab("K") +
  ylab(expression(paste("Calcium ", bgroup("(", frac('mg', '100g'), ")"), sep = ""))) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    text = element_text(size = 6.5))


## Habitat
pd_habitat <- bind_rows(partialPlot(rfm_trip, pred.data = x, x.var = "habitat",
  plot = FALSE))
pd_habitat$x <- fct_reorder(pd_habitat$x, pd_habitat$y, .desc = TRUE)
b <- ggplot(pd_habitat, aes(x = x, y = y)) +
  geom_bar(stat = "identity") +
  # geom_hline(yintercept = dri$dri[dri$nutrients == "calcium_mg"],
  #   linetype = "dashed") +
  coord_cartesian(ylim = c(30, 140)) +
  xlab("Habitat") +
  ylab(expression(paste("Calcium ", bgroup("(", frac('mg', '100g'), ")"), sep = ""))) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    text = element_text(size = 6.5),
    axis.text.x = element_text(angle = 45, hjust = 1))


## Trophic Level
pd_troph <- bind_rows(partialPlot(rfm_trip, pred.data = x, x.var = "troph",
  plot = FALSE, n.pt = 200))
e <- ggplot(pd_troph, aes(x = x, y = y)) +
  geom_line() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  # geom_hline(yintercept = dri$dri[dri$nutrients == "calcium_mg"],
  #   linetype = "dashed") +
  coord_cartesian(ylim = c(30, 140)) +
  xlab("Trophic Level") +
  ylab(expression(paste("Calcium ", bgroup("(", frac('mg', '100g'), ")"), sep = ""))) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    text = element_text(size = 6.5))


## Maturity

# Subset for only L/Lmat < 10
x2 <- filter(x, LLmat < 10)

# Plot results
pd_LLmat <- bind_rows(partialPlot(rfm_trip, pred.data = x2, x.var = "LLmat",
  plot = FALSE, n.pt = 200))
d <- ggplot(pd_LLmat, aes(x = x, y = y)) +
  geom_line() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  # geom_hline(yintercept = dri$dri[dri$nutrients == "calcium_mg"],
  #   linetype = "dashed") +
  coord_cartesian(ylim = c(30, 140)) +
  #xlab(expression(paste("Maturity ", bgroup("(", frac(L, L[mat]), ")")))) +
  xlab("Maturity") +
  ylab(expression(paste("Calcium ", bgroup("(", frac('mg', '100g'), ")"), sep = ""))) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    text = element_text(size = 6.5))


## Length

# Subset for only length < 100 cm
x2 <- filter(x, length_cm < 100)

# Plot results
pd_length <- bind_rows(partialPlot(rfm_trip, pred.data = x2, x.var = "length_cm",
  plot = FALSE, n.pt = 200))
f <- ggplot(pd_length, aes(x = x, y = y)) +
  geom_line() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  # geom_hline(yintercept = dri$dri[dri$nutrients == "calcium_mg"],
  #   linetype = "dashed") +
  coord_cartesian(ylim = c(30, 140)) +
  xlab("Length (cm)") +
  ylab(expression(paste("Calcium ", bgroup("(", frac('mg', '100g'), ")"), sep = ""))) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    text = element_text(size = 6.5))


## Combine plots for figure
ggarrange(a, b, c, d, e, f,
  ncol = 2, nrow = 3,
  labels = "AUTO",
  font.label = list(size = 7))

# Save image
ggsave("figures/06_RandomForestResults_Calcium.pdf",
  device = "pdf",
  width = 172, height = 200,
  units = "mm",
  dpi = 300)




##### RFM 04 - Iron #####

# Select desired variables only
trip2 <- select(trip, c("troph", "k", "habitat", "LLmat", "length_cm", "iron_mg.per.100g"))

# For reproducibility
set.seed(1234)

# Bootstrap sample for model training
trip2.train <- slice_sample(trip2, prop = 1, replace = TRUE)

# The out of bag sample for testing
trip2.test <- setdiff(trip2, trip2.train)

# For reproducibility
set.seed(1234)

# Random forest model with training and testing data
x <- select(trip2.train, -iron_mg.per.100g)
y <- trip2.train$iron_mg.per.100g
xtest <- select(trip2.test, -iron_mg.per.100g)
ytest <- trip2.test$iron_mg.per.100g
rfm_trip <- randomForest(x = x, y = y, xtest = xtest, ytest = ytest,
  mtry = 3, ntree = 501,
  importance = TRUE, proximity = TRUE, keep.forest = TRUE)

# Summary of model
rfm_trip

# 86.26 train / 73.13 test % var explained at mtry = 2
# 86.66 train / 73.17 test % var explained at mtry = 3
# 86.48 train / 72.56 test % var explained at mtry = 4
# 85.96 train / 71.88 test % var explained at mtry = 5

# We chose mtry = 3 because it has the best fit

# Visualize variable importance
ImpData <- as.data.frame(importance(rfm_trip))
ImpData$Var.Names <- row.names(ImpData)
ImpData$Var.Names <- c("Trophic Level", "K", "Habitat",  "Maturity", "Length")
ImpData$Var.Names <- fct_reorder(ImpData$Var.Names, ImpData$`%IncMSE`)
colnames(ImpData)[colnames(ImpData) == "IncNodePurity"] <- "Increase in Node Purity"
a <- ggplot(ImpData, aes(x=Var.Names, y=`%IncMSE`)) +
  geom_segment( aes(x=Var.Names, xend=Var.Names, y=0, yend=`%IncMSE`), color="skyblue") +
  geom_point(aes(size = `Increase in Node Purity`), color="blue", alpha=0.6) +
  theme_light() +
  coord_flip() +
  xlab("") +
  ylab("% Increase in Mean Squared Error") +
  theme_pubr() +
  theme(legend.position = "bottom",
    plot.margin = unit(c(0.5, 1, 0.5, 0), "cm"),
    text = element_text(size = 7),
    axis.text.y = element_text(size = 5),
    axis.text.x = element_text(size = 5))




##### Partial Plots #####

## K
pd_k <- bind_rows(partialPlot(rfm_trip, pred.data = x, x.var = "k",
  plot = FALSE, n.pt = 200))
c <- ggplot(pd_k, aes(x = x, y = y)) +
  geom_line() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  # geom_hline(yintercept = dri$dri[dri$nutrients == "iron_mg"],
  #   linetype = "dashed") +
  coord_cartesian(xlim = c(0, 2), ylim = c(0.5, 2)) +
  xlab("K") +
  ylab(expression(paste("Iron ", bgroup("(", frac('mg', '100g'), ")"), sep = ""))) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    text = element_text(size = 6.5))


## Habitat
pd_habitat <- bind_rows(partialPlot(rfm_trip, pred.data = x, x.var = "habitat",
  plot = FALSE))
pd_habitat$x <- fct_reorder(pd_habitat$x, pd_habitat$y, .desc = TRUE)
b <- ggplot(pd_habitat, aes(x = x, y = y)) +
  geom_bar(stat = "identity") +
  # geom_hline(yintercept = dri$dri[dri$nutrients == "iron_mg"],
  #   linetype = "dashed") +
  coord_cartesian(ylim = c(0.5, 2)) +
  xlab("Habitat") +
  ylab(expression(paste("Iron ", bgroup("(", frac('mg', '100g'), ")"), sep = ""))) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    text = element_text(size = 6.5),
    axis.text.x = element_text(angle = 45, hjust = 1))


## Trophic Level
pd_troph <- bind_rows(partialPlot(rfm_trip, pred.data = x, x.var = "troph",
  plot = FALSE, n.pt = 200))
d <- ggplot(pd_troph, aes(x = x, y = y)) +
  geom_line() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  # geom_hline(yintercept = dri$dri[dri$nutrients == "iron_mg"],
  #   linetype = "dashed") +
  coord_cartesian(ylim = c(0.5, 2)) +
  xlab("Trophic Level") +
  ylab(expression(paste("Iron ", bgroup("(", frac('mg', '100g'), ")"), sep = ""))) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    text = element_text(size = 6.5))


## Maturity

# Subset for only L/Lmat < 10
x2 <- filter(x, LLmat < 10)

# Plot results
pd_LLmat <- bind_rows(partialPlot(rfm_trip, pred.data = x2, x.var = "LLmat",
  plot = FALSE, n.pt = 200))
e <- ggplot(pd_LLmat, aes(x = x, y = y)) +
  geom_line() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  # geom_hline(yintercept = dri$dri[dri$nutrients == "iron_mg"],
  #   linetype = "dashed") +
  coord_cartesian(ylim = c(0.5, 2)) +
  #xlab(expression(paste("Maturity ", bgroup("(", frac(L, L[mat]), ")")))) +
  xlab("Maturity") +
  ylab(expression(paste("Iron ", bgroup("(", frac('mg', '100g'), ")"), sep = ""))) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    text = element_text(size = 6.5))


## Length

# Subset for only length < 100 cm
x2 <- filter(x, length_cm < 100)

# Plot results
pd_length <- bind_rows(partialPlot(rfm_trip, pred.data = x2, x.var = "length_cm",
  plot = FALSE, n.pt = 200))
f <- ggplot(pd_length, aes(x = x, y = y)) +
  geom_line() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  # geom_hline(yintercept = dri$dri[dri$nutrients == "iron_mg"],
  #   linetype = "dashed") +
  coord_cartesian(ylim = c(0.5, 2)) +
  xlab("Length (cm)") +
  ylab(expression(paste("Iron ", bgroup("(", frac('mg', '100g'), ")"), sep = ""))) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    text = element_text(size = 6.5))


## Combine plots for figure
ggarrange(a, b, c, d, e, f,
  ncol = 2, nrow = 3,
  labels = "AUTO",
  font.label = list(size = 7))

# Save image
ggsave("figures/06_RandomForestResults_Iron.pdf",
  device = "pdf",
  width = 172, height = 200,
  units = "mm",
  dpi = 300)




##### RFM 05 - Omega 3 #####

# Select desired variables only
trip2 <- select(trip, c("troph", "k", "habitat", "LLmat", "length_cm", "omega.3_g.per.100g"))

# For reproducibility
set.seed(1234)

# Bootstrap sample for model training
trip2.train <- slice_sample(trip2, prop = 1, replace = TRUE)

# The out of bag sample for testing
trip2.test <- setdiff(trip2, trip2.train)


# Random forest model with training and testing data
x <- select(trip2.train, -omega.3_g.per.100g)
y <- trip2.train$omega.3_g.per.100g
xtest <- select(trip2.test, -omega.3_g.per.100g)
ytest <- trip2.test$omega.3_g.per.100g
set.seed(1234)
rfm_trip <- randomForest(x = x, y = y, xtest = xtest, ytest = ytest,
  mtry = 3, ntree = 501,
  importance = TRUE, proximity = TRUE, keep.forest = TRUE)

# Summary of model
rfm_trip

# 88.93 train / 80.57 test % var explained at mtry = 2
# 89.49 train / 80.31 test % var explained at mtry = 3
# 89.44 train / 79.43 test % var explained at mtry = 4
# 89.63 / 78.17 test % var explained at mtry = 5

# We chose mtry = 3 because it has the best fit for testing and training data

# Visualize variable importance
ImpData <- as.data.frame(importance(rfm_trip))
ImpData$Var.Names <- row.names(ImpData)
ImpData$Var.Names <- c("Trophic Level", "K", "Habitat",  "Maturity", "Length")
ImpData$Var.Names <- fct_reorder(ImpData$Var.Names, ImpData$`%IncMSE`)
colnames(ImpData)[colnames(ImpData) == "IncNodePurity"] <- "Increase in Node Purity"
a <- ggplot(ImpData, aes(x=Var.Names, y=`%IncMSE`)) +
  geom_segment( aes(x=Var.Names, xend=Var.Names, y=0, yend=`%IncMSE`), color="skyblue") +
  geom_point(aes(size = `Increase in Node Purity`), color="blue", alpha=0.6) +
  theme_light() +
  coord_flip() +
  xlab("") +
  ylab("% Increase in Mean Squared Error") +
  theme_pubr() +
  theme(legend.position = "bottom",
    text = element_text(size = 7),
    axis.text.y = element_text(size = 5),
    axis.text.x = element_text(size = 5))




##### Partial Plots #####

## K
pd_k <- bind_rows(partialPlot(rfm_trip, pred.data = x, x.var = "k",
  plot = FALSE, n.pt = 200))
d <- ggplot(pd_k, aes(x = x, y = y)) +
  geom_line() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  # geom_hline(yintercept = dri$dri[dri$nutrients == "omega3_g"],
  #   linetype = "dashed") +
  coord_cartesian(xlim = c(0, 2), ylim = c(0.1, 0.3)) +
  xlab("K") +
  ylab(expression(paste("Omega-3 ", bgroup("(", frac('g', '100g'), ")"), sep = ""))) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    text = element_text(size = 6.5))


## Habitat
pd_habitat <- bind_rows(partialPlot(rfm_trip, pred.data = x, x.var = "habitat",
  plot = FALSE))
pd_habitat$x <- fct_reorder(pd_habitat$x, pd_habitat$y, .desc = TRUE)
b <- ggplot(pd_habitat, aes(x = x, y = y)) +
  geom_bar(stat = "identity") +
  # geom_hline(yintercept = dri$dri[dri$nutrients == "omega3_g"],
  #   linetype = "dashed") +
  coord_cartesian(ylim = c(0.1, 0.3)) +
  xlab("Habitat") +
  ylab(expression(paste("Omega-3 ", bgroup("(", frac('g', '100g'), ")"), sep = ""))) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    text = element_text(size = 6.5),
    axis.text.x = element_text(angle = 45, hjust = 1))


## Trophic Level
pd_troph <- bind_rows(partialPlot(rfm_trip, pred.data = x, x.var = "troph",
  plot = FALSE, n.pt = 200))
c <- ggplot(pd_troph, aes(x = x, y = y)) +
  geom_line() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  # geom_hline(yintercept = dri$dri[dri$nutrients == "omega3_g"],
  #   linetype = "dashed") +
  coord_cartesian(ylim = c(0.1, 0.3)) +
  xlab("Trophic Level") +
  ylab(expression(paste("Omega-3 ", bgroup("(", frac('g', '100g'), ")"), sep = ""))) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    text = element_text(size = 6.5))


## Maturity

# Subset for only L/Lmat < 10
x2 <- filter(x, LLmat < 10)

# Plot results
pd_LLmat <- bind_rows(partialPlot(rfm_trip, pred.data = x2, x.var = "LLmat",
  plot = FALSE, n.pt = 200))
e <- ggplot(pd_LLmat, aes(x = x, y = y)) +
  geom_line() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  # geom_hline(yintercept = dri$dri[dri$nutrients == "omega3_g"],
  #   linetype = "dashed") +
  coord_cartesian(ylim = c(0.1, 0.3)) +
  #xlab(expression(paste("Maturity ", bgroup("(", frac(L, L[mat]), ")")))) +
  xlab("Maturity") +
  ylab(expression(paste("Omega-3 ", bgroup("(", frac('g', '100g'), ")"), sep = ""))) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    text = element_text(size = 6.5))


## Length

# Subset for only length < 100 cm
x2 <- filter(x, length_cm < 100)

# Plot results
pd_length <- bind_rows(partialPlot(rfm_trip, pred.data = x2, x.var = "length_cm",
  plot = FALSE, n.pt = 200))
f <- ggplot(pd_length, aes(x = x, y = y)) +
  geom_line() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  # geom_hline(yintercept = dri$dri[dri$nutrients == "omega3_g"],
  #   linetype = "dashed") +
  coord_cartesian(ylim = c(0.1, 0.3)) +
  xlab("Length (cm)") +
  ylab(expression(paste("Omega-3 ", bgroup("(", frac('g', '100g'), ")"), sep = ""))) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    text = element_text(size = 6.5))


## Combine plots for figure
ggarrange(a, b, c, d, e, f,
  ncol = 2, nrow = 3,
  labels = "AUTO",
  font.label = list(size = 7))

# Save image
ggsave("figures/06_RandomForestResults_Omega3.pdf",
  device = "pdf",
  width = 172, height = 200,
  units = "mm",
  dpi = 300)




##### RFM 06 - Selenium #####

# Select desired variables only
trip2 <- select(trip, c("troph", "k", "habitat", "LLmat", "length_cm", "selenium_ug.per.100g"))

# For reproducibility
set.seed(1234)

# Bootstrap sample for model training
trip2.train <- slice_sample(trip2, prop = 1, replace = TRUE)

# The out of bag sample for testing
trip2.test <- setdiff(trip2, trip2.train)

# Random forest model with training and testing data
x <- select(trip2.train, -selenium_ug.per.100g)
y <- trip2.train$selenium_ug.per.100g
xtest <- select(trip2.test, -selenium_ug.per.100g)
ytest <- trip2.test$selenium_ug.per.100g
set.seed(1234)
rfm_trip <- randomForest(x = x, y = y, xtest = xtest, ytest = ytest,
  mtry = 3, ntree = 501,
  importance = TRUE, proximity = TRUE, keep.forest = TRUE)

# Summary of model
rfm_trip

# 89.11 train / 72.87 test % var explained at mtry = 2
# 89.23 train / 72.44 test % var explained at mtry = 3
# 89.32 train / 72.07 test % var explained at mtry = 4
# 89.14 train / 70.77 test % var explained at mtry = 5

# We chose mtry = 3 because it has the best fit for training and testing data

# Visualize variable importance
ImpData <- as.data.frame(importance(rfm_trip))
ImpData$Var.Names <- row.names(ImpData)
ImpData$Var.Names <- c("Trophic Level", "K", "Habitat",  "Maturity", "Length")
ImpData$Var.Names <- fct_reorder(ImpData$Var.Names, ImpData$`%IncMSE`)
colnames(ImpData)[colnames(ImpData) == "IncNodePurity"] <- "Increase in Node Purity"
a <- ggplot(ImpData, aes(x=Var.Names, y=`%IncMSE`)) +
  geom_segment( aes(x=Var.Names, xend=Var.Names, y=0, yend=`%IncMSE`), color="skyblue") +
  geom_point(aes(size = `Increase in Node Purity`), color="blue", alpha=0.6) +
  theme_light() +
  coord_flip() +
  xlab("") +
  ylab("% Increase in Mean Squared Error") +
  theme_pubr() +
  theme(legend.position = "bottom",
    text = element_text(size = 7),
    axis.text.y = element_text(size = 5),
    axis.text.x = element_text(size = 5))




##### Partial Plots #####

## K
pd_k <- bind_rows(partialPlot(rfm_trip, pred.data = x, x.var = "k",
  plot = FALSE, n.pt = 200))
b <- ggplot(pd_k, aes(x = x, y = y)) +
  geom_line() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  geom_hline(yintercept = dri$dri[dri$nutrients == "selenium_ug"],
    linetype = "dashed") +
  coord_cartesian(xlim = c(0, 2), ylim = c(0, 150)) +
  xlab("K") +
  ylab(expression(paste("Selenium ", bgroup("(", frac(paste("\u00b5", g, sep = ""), '100g'), ")"), sep = ""))) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    text = element_text(size = 6.5))


## Habitat
pd_habitat <- bind_rows(partialPlot(rfm_trip, pred.data = x, x.var = "habitat",
  plot = FALSE))
pd_habitat$x <- fct_reorder(pd_habitat$x, pd_habitat$y, .desc = TRUE)
e <- ggplot(pd_habitat, aes(x = x, y = y)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = dri$dri[dri$nutrients == "selenium_ug"],
    linetype = "dashed") +
  coord_cartesian(ylim = c(0, 150)) +
  xlab("Habitat") +
  ylab(expression(paste("Selenium ", bgroup("(", frac(paste("\u00b5", g, sep = ""), '100g'), ")"), sep = ""))) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    text = element_text(size = 6.5),
    axis.text.x = element_text(angle = 45, hjust = 1))


## Trophic Level
pd_troph <- bind_rows(partialPlot(rfm_trip, pred.data = x, x.var = "troph",
  plot = FALSE, n.pt = 200))
f <- ggplot(pd_troph, aes(x = x, y = y)) +
  geom_line() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  geom_hline(yintercept = dri$dri[dri$nutrients == "selenium_ug"],
    linetype = "dashed") +
  coord_cartesian(ylim = c(0, 150)) +
  xlab("Trophic Level") +
  ylab(expression(paste("Selenium ", bgroup("(", frac(paste("\u00b5", g, sep = ""), '100g'), ")"), sep = ""))) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.9, 0.5), "cm"),
    text = element_text(size = 6.5),
    axis.title.x = element_text(margin = unit(c(0.8, 0, 0, 0), "cm")))


## Maturity

# Subset for only L/Lmat < 10
x2 <- filter(x, LLmat < 10)

# Plot results
pd_LLmat <- bind_rows(partialPlot(rfm_trip, pred.data = x2, x.var = "LLmat",
  plot = FALSE, n.pt = 200))
c <- ggplot(pd_LLmat, aes(x = x, y = y)) +
  geom_line() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  geom_hline(yintercept = dri$dri[dri$nutrients == "selenium_ug"],
    linetype = "dashed") +
  coord_cartesian(ylim = c(0, 150)) +
  #xlab(expression(paste("Maturity ", bgroup("(", frac(L, L[mat]), ")")))) +
  xlab("Maturity") +
  ylab(expression(paste("Selenium ", bgroup("(", frac(paste("\u00b5", g, sep = ""), '100g'), ")"), sep = ""))) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    text = element_text(size = 6.5))


## Length

# Subset for only length < 100 cm
x2 <- filter(x, length_cm < 100)

# Plot results
pd_length <- bind_rows(partialPlot(rfm_trip, pred.data = x2, x.var = "length_cm",
  plot = FALSE, n.pt = 200))
d <- ggplot(pd_length, aes(x = x, y = y)) +
  geom_line() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  geom_hline(yintercept = dri$dri[dri$nutrients == "selenium_ug"],
    linetype = "dashed") +
  coord_cartesian(ylim = c(0, 150)) +
  xlab("Length (cm)") +
  ylab(expression(paste("Selenium ", bgroup("(", frac(paste("\u00b5", g, sep = ""), '100g'), ")"), sep = ""))) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    text = element_text(size = 6.5))


## Combine plots for figure
ggarrange(a, b, c, d, e, f,
  ncol = 2, nrow = 3,
  labels = "AUTO",
  font.label = list(size = 7))

# Save image
ggsave("figures/06_RandomForestResults_Selenium.pdf",
  device = "pdf",
  width = 172, height = 200,
  units = "mm",
  dpi = 300)




##### RFM 07 - Vitamin A #####

# Select desired variables only
trip2 <- select(trip, c("troph", "k", "habitat", "LLmat", "length_cm", "vitamin.a_ug.per.100g"))

# For reproducibility
set.seed(1234)

# Bootstrap sample for model training
trip2.train <- slice_sample(trip2, prop = 1, replace = TRUE)

# The out of bag sample for testing
trip2.test <- setdiff(trip2, trip2.train)

# Random forest model with training and testing data
x <- select(trip2.train, -vitamin.a_ug.per.100g)
y <- trip2.train$vitamin.a_ug.per.100g
xtest <- select(trip2.test, -vitamin.a_ug.per.100g)
ytest <- trip2.test$vitamin.a_ug.per.100g
set.seed(1234)
rfm_trip <- randomForest(x = x, y = y, xtest = xtest, ytest = ytest,
  mtry = 4, ntree = 501,
  importance = TRUE, proximity = TRUE, keep.forest = TRUE)

# Summary of model
rfm_trip

# 82.78 train / 59.48 test % var explained at mtry = 2
# 83.71 train / 60.28 test % var explained at mtry = 3
# 83.82 train / 60.43 test % var explained at mtry = 4
# 83.88 train / 59.41 test % var explained at mtry = 5

# We chose mtry = 4 because it has the best fit for training and testing data

# Visualize variable importance
ImpData <- as.data.frame(importance(rfm_trip))
ImpData$Var.Names <- row.names(ImpData)
ImpData$Var.Names <- c("Trophic Level", "K", "Habitat",  "Maturity", "Length")
ImpData$Var.Names <- fct_reorder(ImpData$Var.Names, ImpData$`%IncMSE`)
colnames(ImpData)[colnames(ImpData) == "IncNodePurity"] <- "Increase in Node Purity"
a <- ggplot(ImpData, aes(x=Var.Names, y=`%IncMSE`)) +
  geom_segment( aes(x=Var.Names, xend=Var.Names, y=0, yend=`%IncMSE`), color="skyblue") +
  geom_point(aes(size = `Increase in Node Purity`), color="blue", alpha=0.6) +
  theme_light() +
  coord_flip() +
  xlab("") +
  ylab("% Increase in Mean Squared Error\n\n\nIncrease in Node Purity") +
  theme_pubr() +
  theme(legend.position = c(.5, -0.32), legend.direction = "horizontal",
    legend.title = element_blank(),
    plot.margin = unit(c(0.5, 0.5, 1, 0.5), "cm"),
    text = element_text(size = 7),
    axis.text.y = element_text(size = 5),
    axis.text.x = element_text(size = 5))




##### Partial Plots #####

## K
pd_k <- bind_rows(partialPlot(rfm_trip, pred.data = x, x.var = "k",
  plot = FALSE, n.pt = 200))
e <- ggplot(pd_k, aes(x = x, y = y)) +
  geom_line() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  # geom_hline(yintercept = dri$dri[dri$nutrients == "vitamina_ug"],
  #   linetype = "dashed") +
  coord_cartesian(xlim = c(0, 2), ylim = c(0, 100)) +
  xlab("K") +
  ylab(expression(paste("Vitamin A ", bgroup("(", frac(paste("\u00b5", g, sep = ""), '100g'), ")"), sep = ""))) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    text = element_text(size = 6.5))


## Habitat
pd_habitat <- bind_rows(partialPlot(rfm_trip, pred.data = x, x.var = "habitat",
  plot = FALSE))
pd_habitat$x <- fct_reorder(pd_habitat$x, pd_habitat$y, .desc = TRUE)
c <- ggplot(pd_habitat, aes(x = x, y = y)) +
  geom_bar(stat = "identity") +
  # geom_hline(yintercept = dri$dri[dri$nutrients == "vitamina_ug"],
  #   linetype = "dashed") +
  coord_cartesian(ylim = c(0, 100)) +
  xlab("Habitat") +
  ylab(expression(paste("Vitamin A ", bgroup("(", frac(paste("\u00b5", g, sep = ""), '100g'), ")"), sep = ""))) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    text = element_text(size = 6.5),
    axis.text.x = element_text(angle = 45, hjust = 1))


## Trophic Level
pd_troph <- bind_rows(partialPlot(rfm_trip, pred.data = x, x.var = "troph",
  plot = FALSE, n.pt = 200))
b <- ggplot(pd_troph, aes(x = x, y = y)) +
  geom_line() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  # geom_hline(yintercept = dri$dri[dri$nutrients == "vitamina_ug"],
  #   linetype = "dashed") +
  coord_cartesian(ylim = c(0, 100)) +
  xlab("Trophic Level") +
  ylab(expression(paste("Vitamin A ", bgroup("(", frac(paste("\u00b5", g, sep = ""), '100g'), ")"), sep = ""))) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    text = element_text(size = 6.5))


## Maturity

# Subset for only L/Lmat < 10
x2 <- filter(x, LLmat < 10)

# Plot results
pd_LLmat <- bind_rows(partialPlot(rfm_trip, pred.data = x2, x.var = "LLmat",
  plot = FALSE, n.pt = 200))
f <- ggplot(pd_LLmat, aes(x = x, y = y)) +
  geom_line() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  # geom_hline(yintercept = dri$dri[dri$nutrients == "vitamina_ug"],
  #   linetype = "dashed") +
  coord_cartesian(ylim = c(0, 100)) +
  #xlab(expression(paste("Maturity ", bgroup("(", frac(L, L[mat]), ")")))) +
  xlab("Maturity") +
  ylab(expression(paste("Vitamin A ", bgroup("(", frac(paste("\u00b5", g, sep = ""), '100g'), ")"), sep = ""))) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    text = element_text(size = 6.5))


## Length

# Subset for only length < 100 cm
x2 <- filter(x, length_cm < 100)

# Plot results
pd_length <- bind_rows(partialPlot(rfm_trip, pred.data = x2, x.var = "length_cm",
  plot = FALSE, n.pt = 200))
d <- ggplot(pd_length, aes(x = x, y = y)) +
  geom_line() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  # geom_hline(yintercept = dri$dri[dri$nutrients == "vitamina_ug"],
  #   linetype = "dashed") +
  coord_cartesian(ylim = c(0, 100)) +
  xlab("Length (cm)") +
  ylab(expression(paste("Vitamin A ", bgroup("(", frac(paste("\u00b5", g, sep = ""), '100g'), ")"), sep = ""))) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.9, 0.5), "cm"),
    text = element_text(size = 6.5),
    axis.title.x = element_text(margin = unit(c(0.8, 0, 0, 0), "cm")))


## Combine plots for figure
ggarrange(a, b, c, d, e, f,
  ncol = 2, nrow = 3,
  labels = "AUTO",
  font.label = list(size = 7))

# Save image
ggsave("figures/06_RandomForestResults_VitaminA.pdf",
  device = "pdf",
  width = 172, height = 200,
  units = "mm",
  dpi = 300)





##### RFM 08 - Zinc #####

# Select desired variables only
trip2 <- select(trip, c("troph", "k", "habitat", "LLmat", "length_cm", "zinc_mg.per.100g"))

# For reproducibility
set.seed(1234)

# Bootstrap sample for model training
trip2.train <- slice_sample(trip2, prop = 1, replace = TRUE)

# The out of bag sample for testing
trip2.test <- setdiff(trip2, trip2.train)

# Random forest model with training and testing data
x <- select(trip2.train, -zinc_mg.per.100g)
y <- trip2.train$zinc_mg.per.100g
xtest <- select(trip2.test, -zinc_mg.per.100g)
ytest <- trip2.test$zinc_mg.per.100g
set.seed(1234)
rfm_trip <- randomForest(x = x, y = y, xtest = xtest, ytest = ytest,
  mtry = 4, ntree = 501,
  importance = TRUE, proximity = TRUE, keep.forest = TRUE)

# Summary of model
rfm_trip

# 85.5 train / 66.52 test % var explained at mtry = 2
# 86.28 train / 67.28 test % var explained at mtry = 3
# 86.54 train / 67.53 test % var explained at mtry = 4
# 86.5 train / 67.28 test % var explained at mtry = 5

# We chose mtry = 4 because it has the best fit for training and testing data

# Visualize variable importance
ImpData <- as.data.frame(importance(rfm_trip))
ImpData$Var.Names <- row.names(ImpData)
ImpData$Var.Names <- c("Trophic Level", "K", "Habitat",  "Maturity", "Length")
ImpData$Var.Names <- fct_reorder(ImpData$Var.Names, ImpData$`%IncMSE`)
colnames(ImpData)[colnames(ImpData) == "IncNodePurity"] <- "Increase in Node Purity"
a <- ggplot(ImpData, aes(x=Var.Names, y=`%IncMSE`)) +
  geom_segment( aes(x=Var.Names, xend=Var.Names, y=0, yend=`%IncMSE`), color="skyblue") +
  geom_point(aes(size = `Increase in Node Purity`), color="blue", alpha=0.6) +
  theme_light() +
  coord_flip() +
  xlab("") +
  ylab("% Increase in Mean Squared Error\n\nIncrease in Node Purity") +
  theme_pubr() +
  theme(legend.position = c(.5, -0.45), legend.direction = "horizontal",
    legend.title = element_blank(),
    plot.margin = unit(c(0.5, 0.5, 1, 0.5), "cm"),
    text = element_text(size = 7),
    axis.text.y = element_text(size = 5),
    axis.text.x = element_text(size = 5))



##### Partial Plots #####

## K
pd_k <- bind_rows(partialPlot(rfm_trip, pred.data = x, x.var = "k",
  plot = FALSE, n.pt = 200))
c <- ggplot(pd_k, aes(x = x, y = y)) +
  geom_line() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  geom_hline(yintercept = dri$dri[dri$nutrients == "zinc_mg"],
    linetype = "dashed") +
  coord_cartesian(xlim = c(0, 2), ylim = c(0, 3.1)) +
  xlab("K") +
  ylab(expression(paste("Zinc ", bgroup("(", frac('mg', '100g'), ")"), sep = ""))) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    text = element_text(size = 6.5))


## Habitat
pd_habitat <- bind_rows(partialPlot(rfm_trip, pred.data = x, x.var = "habitat",
  plot = FALSE))
pd_habitat$x <- fct_reorder(pd_habitat$x, pd_habitat$y, .desc = TRUE)
f <- ggplot(pd_habitat, aes(x = x, y = y)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = dri$dri[dri$nutrients == "zinc_mg"],
    linetype = "dashed") +
  coord_cartesian(ylim = c(0, 3.1)) +
  xlab("Habitat") +
  ylab(expression(paste("Zinc ", bgroup("(", frac('mg', '100g'), ")"), sep = ""))) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    text = element_text(size = 6.5),
    axis.text.x = element_text(angle = 45, hjust = 1))


## Trophic Level
pd_troph <- bind_rows(partialPlot(rfm_trip, pred.data = x, x.var = "troph",
  plot = FALSE, n.pt = 200))
b <- ggplot(pd_troph, aes(x = x, y = y)) +
  geom_line() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  geom_hline(yintercept = dri$dri[dri$nutrients == "zinc_mg"],
    linetype = "dashed") +
  coord_cartesian(ylim = c(0, 3.1)) +
  xlab("Trophic Level") +
  ylab(expression(paste("Zinc ", bgroup("(", frac('mg', '100g'), ")"), sep = ""))) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    text = element_text(size = 6.5))


## Maturity

# Subset for only L/Lmat < 10
x2 <- filter(x, LLmat < 10)

# Plot results
pd_LLmat <- bind_rows(partialPlot(rfm_trip, pred.data = x2, x.var = "LLmat",
  plot = FALSE, n.pt = 200))
e <- ggplot(pd_LLmat, aes(x = x, y = y)) +
  geom_line() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  geom_hline(yintercept = dri$dri[dri$nutrients == "zinc_mg"],
    linetype = "dashed") +
  coord_cartesian(ylim = c(0, 3.1)) +
  #xlab(expression(paste("Maturity ", bgroup("(", frac(L, L[mat]), ")")))) +
  xlab("Maturity") +
  ylab(expression(paste("Zinc ", bgroup("(", frac('mg', '100g'), ")"), sep = ""))) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.9, 0.5), "cm"),
    text = element_text(size = 6.5),
    axis.title.x = element_text(margin = unit(c(0.8, 0, 0, 0), "cm")))


## Length

# Subset for only length < 100 cm
x2 <- filter(x, length_cm < 100)

# Plot results
pd_length <- bind_rows(partialPlot(rfm_trip, pred.data = x2, x.var = "length_cm",
  plot = FALSE, n.pt = 200))
d <- ggplot(pd_length, aes(x = x, y = y)) +
  geom_line() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  geom_hline(yintercept = dri$dri[dri$nutrients == "zinc_mg"],
    linetype = "dashed") +
  coord_cartesian(ylim = c(0, 3.1)) +
  xlab("Length (cm)") +
  ylab(expression(paste("Zinc ", bgroup("(", frac('mg', '100g'), ")"), sep = ""))) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    text = element_text(size = 6.5))


## Combine plots for figure
ggarrange(a, b, c, d, e, f,
  ncol = 2, nrow = 3,
  labels = "AUTO",
  font.label = list(size = 7))

# Save image
ggsave("figures/06_RandomForestResults_Zinc.pdf",
  device = "pdf",
  width = 172, height = 200,
  units = "mm",
  dpi = 300)



