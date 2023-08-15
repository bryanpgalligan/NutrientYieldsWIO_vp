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

# Get variable importance
ImpData <- as.data.frame(importance(rfm_trip))

# Visualize variable importance
ImpData$Var.Names <- row.names(ImpData)
ImpData$Var.Names <- c("Country", "Site", "Gear", "Landings", "CPUE", "Maturity",
  "Length : Optimum Length", "Length", "Species Richness", "Management", "Current",
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
ggsave("figures/06_RandomForestAllVariables.pdf",
  device = "pdf",
  width = 180,
  units = "mm",
  dpi = 300)

# Retain only those variables that performed better than random numbers
trip <- select(trip, c(trip, date, latitude, longitude, source, nutrients_pdv, species_richness, cpue_kg.effort.day, LLmat, gear, LLopt, length_cm))


## Spearman correlation matrix for numeric variables

# Subset of numeric data
x <- select(trip, species_richness:length_cm)
x <- select(x, -gear)

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

# Alternative model using LLmat instead of LLopt
trip.alt <- select(trip, -LLopt)

# L/Lmat ~ L/Lopt >> remove L/Lmat for the main model
trip <- select(trip, -LLmat)




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
  mtry = 2, ntree = 501,
  importance = TRUE, proximity = TRUE, keep.forest = TRUE)

# Summary of model
rfm_trip

# Visualize variable importance
ImpData <- as.data.frame(importance(rfm_trip))
ImpData$Var.Names <- row.names(ImpData)
ImpData$Var.Names <- c("Species Richness", "CPUE", "Gear", "Length : Optimum Length", "Length")
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
    text = element_text(size = 6.5),
    axis.title.x = element_text(vjust = -3))

# # Save plot
# ggsave("figures/06_RandomForestSupp2.pdf",
#   device = "pdf",
#   width = 180,
#   units = "mm",
#   dpi = 300)




##### Partial Plots #####

# Plot the predicted relationship between LLopt and total nutrients
pd_LLopt <- bind_rows(partialPlot(rfm_trip, pred.data = x, x.var = "LLopt",
  plot = FALSE, n.pt = 200))
b <- ggplot(pd_LLopt, aes(x = x, y = 100*y)) +
  geom_line() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  coord_cartesian(xlim = c(0, 5), ylim = c(190, 260)) +
  xlab(expression(paste("Length : Optimum Length ", bgroup("(", frac(L, L[opt]), ")")))) +
  ylab("Nutrient Density (% DRI)") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    text = element_text(size = 6.5))

# Plot the predicted relationship between gear and total nutrients
pd_gear <- bind_rows(partialPlot(rfm_trip, pred.data = x, x.var = "gear",
  plot = FALSE))
pd_gear$x <- fct_reorder(pd_gear$x, pd_gear$y, .desc = TRUE)
pd_gear <- filter(pd_gear, x != "other")
c <- ggplot(pd_gear, aes(x = x, y = 100*y)) +
  geom_bar(stat = "identity") +
  coord_cartesian(ylim = c(190, 260)) +
  xlab("Gear") +
  ylab("Nutrient Density (% DRI)") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    text = element_text(size = 6.5),
    axis.text.x = element_text(angle = 45, hjust = 1))

# Plot the predicted relationship between length and total nutrients
pd_length <- bind_rows(partialPlot(rfm_trip, pred.data = x, x.var = "length_cm",
  plot = FALSE, n.pt = 200))
d <- ggplot(pd_length, aes(x = x, y = 100*y)) +
  geom_line() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  coord_cartesian(xlim = c(0, 100), ylim = c(190, 260)) +
  xlab("Length (cm)") +
  ylab("Nutrient Density (% DRI)") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(plot.margin = unit(c(5, 5, 5, 5), "mm"),
    text = element_text(size = 6.5),
    axis.title.x = element_text(margin = unit(c(10, 0, 0, 0), units = "mm")))

# Plot the predicted relationship between species richness and total nutrients
pd_spr <- bind_rows(partialPlot(rfm_trip, pred.data = x, x.var = "species_richness",
  plot = FALSE, n.pt = 200))
e <- ggplot(pd_spr, aes(x = x, y = 100*y)) +
  geom_line() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  coord_cartesian(xlim = c(0, 5), ylim = c(190, 260)) +
  xlab("Species Richness (per kg)") +
  ylab("Nutrient Density (% DRI)") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    text = element_text(size = 6.5))

# Plot the predicted relationship between CPUE and total nutrients
pd_cpue <- bind_rows(partialPlot(rfm_trip, pred.data = x, x.var = "cpue_kg.effort.day",
  plot = FALSE, n.pt = 200))
f <- ggplot(pd_cpue, aes(x = x, y = 100*y)) +
  geom_line() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  coord_cartesian(xlim = c(0, 10), ylim = c(190, 260)) +
  xlab(expression(paste("Catch per Unit Effort ", bgroup("(", frac(kg, "fisher x day"), ")")))) +
  ylab("Nutrient Density (% DRI)") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    text = element_text(size = 6.5))

# Combine plots for fig 1
ggarrange(a, b, c, d,
  ncol = 2, nrow = 2,
  labels = "AUTO",
  font.label = list(size = 7))

# Save image
ggsave("figures/06_RandomForestResults.pdf",
  device = "pdf",
  width = 174, height = 120,
  units = "mm",
  dpi = 300)

# Combine remaining partial plots for supp mat
ggarrange(e, f, ncol = 1, nrow = 2, labels = "AUTO", font.label = list(size = 7))

# Save image
ggsave("figures/06_RandomForestAdditionalResults.pdf",
  device = "pdf",
  width = 90,
  units = "mm",
  dpi = 300)




##### RFM 02 - Pruned #####

## Evaluate partial dependence plots for shape-based significance (Greenwell et al. 2018).

# L/Lopt = 0.06537712
sd(pd_LLopt$y)

# Gear = 0.1225924
(max(pd_gear$y) - min(pd_gear$y)) / 4

# Length = 0.1121959
sd(pd_length$y)

# Species richness = 0.007335326
sd(pd_spr$y)

# CPUE = 0.00563308
sd(pd_cpue$y)

# Species richness and CPUE are both below 0.01 so are very flat


## A new RFM with species richness and CPUE removed

# Remove species richness and CPUE
trip2 <- select(trip2, -c("species_richness", "cpue_kg.effort.day"))

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
  mtry = 2, ntree = 501,
  importance = TRUE, proximity = TRUE, keep.forest = TRUE)

# Summary of model
rfm_trip

# The new model explains 64% of training data variance and 21% of test data variance

# Visualize variable importance
ImpData <- as.data.frame(importance(rfm_trip))
ImpData$Var.Names <- row.names(ImpData)
#ImpData$Var.Names <- c("Species Richness", "CPUE", "Gear", "Length : Optimum Length", "Length")
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
    text = element_text(size = 7),
    axis.text.y = element_text(size = 5),
    axis.text.x = element_text(size = 5))




##### Alternate (L/Lmat) Model #####

# Data for analysis
trip3 <- select(trip.alt, -c("trip", "date", "latitude", "longitude", "source"))

# For reproducibility
set.seed(1234)

# Impute NAs using random forest proximity
trip3 <- rfImpute(nutrients_pdv ~ ., data = trip3)

# For reproducibility
set.seed(1234)

# Bootstrap sample for model training
trip3.train <- slice_sample(trip3, prop = 1, replace = TRUE)

# The out of bag sample for testing
trip3.test <- setdiff(trip3, trip3.train)

# For reproducibility
set.seed(1234)

# Random forest model with training and testing data
x <- select(trip3.train, -nutrients_pdv)
y <- trip3.train$nutrients_pdv
xtest <- select(trip3.test, -nutrients_pdv)
ytest <- trip3.test$nutrients_pdv
rfm_trip2 <- randomForest(x = x, y = y, xtest = xtest, ytest = ytest,
  mtry = 2, ntree = 501,
  importance = TRUE, proximity = TRUE, keep.forest = TRUE)

# Summary of model
rfm_trip2

# Visualize variable importance
ImpData <- as.data.frame(importance(rfm_trip2))
ImpData$Var.Names <- row.names(ImpData)
ImpData$Var.Names <- c("Species Richness", "CPUE", "Maturity", "Gear", "Length")
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
    text = element_text(size = 6.5))

# Plot the predicted relationship between LLmat and total nutrients
pd_LLmat <- bind_rows(partialPlot(rfm_trip2, pred.data = x, x.var = "LLmat",
  plot = FALSE, n.pt = 200))
b <- ggplot(pd_LLmat, aes(x = x, y = 100*y)) +
  geom_line() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  coord_cartesian(xlim = c(0, 5), ylim = c(190, 260)) +
  xlab(expression(paste("Maturity ", bgroup("(", frac(L, L[mat]), ")")))) +
  ylab("Nutrient Density (% DRI)") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    text = element_text(size = 6.5))

# Combine plots for supplementary material
ggarrange(a, b,
  labels = "AUTO",
  font.label = list(size = 7))

# Save image
ggsave("figures/06_RandomForestMaturityResults.pdf",
  device = "pdf",
  width = 180,
  height = 60,
  units = "mm",
  dpi = 300)
