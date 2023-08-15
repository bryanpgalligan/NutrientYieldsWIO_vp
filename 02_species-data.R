## This script retrieves, cleans, and adds to the fish nutrients data from FishBase.

# Retrieve FishBase "estimates," "species," and "taxa" tables
species_data <- left_join(estimate(), species(), by = "SpecCode")
species_data<- left_join(species_data, load_taxa(), by = "SpecCode")

# Select only desired variables
species_data <- select(species_data, Species, Class, Order, Family, Genus.y, MaxLengthTL, Troph, a, b, K, TempPrefMin, TempPrefMean,
  TempPrefMax, FeedingPath, FeedingPathLevel, Calcium, Iron, Omega3, Selenium, VitaminA, Zinc)

# Rename variables
species_data <- rename(species_data,
  c(species = Species,
    class = Class,
    order = Order,
    family = Family,
    genus = Genus.y,
    max.length.tl = MaxLengthTL,
    troph = Troph,
    temp.pref.min = TempPrefMin,
    temp.pref.mean = TempPrefMean,
    temp.pref.max = TempPrefMax,
    feeding.path = FeedingPath,
    feeding.path.level = FeedingPathLevel,
    calcium_mg.100g = Calcium,
    iron_mg.100g = Iron,
    omega3_g.100g = Omega3,
    selenium_ug.100g = Selenium,
    vitamin.a_ug.100g = VitaminA,
    zinc_mg.100g = Zinc
  ))


# Dietary reference intakes (DRI) (recommended dietary allowances or adequate intakes)
# of key nutrients for children 1-3 years old (IOM 2006, 2011)

# List of nutrients
nutrients <- c(
  "calcium_mg",
  "iron_mg",
  "omega3_g",
  "protein_g",
  "selenium_ug",
  "vitamina_ug",
  "zinc_mg"
  )

# List of DRIs for children 1-3 years extracted from IOM (2006, 2011)
dri <- c(
  700, # calcium recommended daily allowance (RDA)
  7, # iron recommended daily allowance (RDA)
  0.7, # alpha-linolenic acid (n-3) adequate intake (AI)
  13, # protein recommended daily allowance (RDA)
  20, # selenium recommended daily allowance (RDA)
  300, # vitamin A recommended daily allowance (RDA)
  3 # zinc recommended daily allowance (RDA)
)

# Make a data frame
dri_nutrients <- data.frame(nutrients, dri)

# Empty column in species data
species_data$total.nutrients <- NA

# Add overall nutrient density scores to species_data
for (i in 1:nrow(species_data)){
  
  if (is.na(species_data$calcium_mg.100g[i]) == FALSE){
  
    # Calculate nutrients per daily value
    ca <- species_data$calcium_mg.100g[i] / dri_nutrients$dri[1]
    fe <- species_data$iron_mg.100g[i] / dri_nutrients$dri[2]
    o3 <- species_data$omega3_g.100g[i] / dri_nutrients$dri[3]
    se <- species_data$selenium_ug.100g[i] / dri_nutrients$dri[5]
    va <- species_data$vitamin.a_ug.100g[i] / dri_nutrients$dri[6]
    zn <- species_data$zinc_mg.100g[i] / dri_nutrients$dri[7]
    
    # Cap values at 1
    if (ca > 1) {ca <- 1}
    if (fe > 1) {fe <- 1}
    if (o3 > 1) {o3 <- 1}
    if (se > 1) {se <- 1}
    if (va > 1) {va <- 1}
    if (zn > 1) {zn <- 1}
    
    # Save to species_data
    species_data$total.nutrients[i] <- ca + fe + o3 + se + va + zn
  }
  
}


## Save files

# Save species data
write.csv(species_data, file = "data/temp-data/02_SpeciesData.csv",
  row.names = FALSE)

# Save dietary reference intakes
write.csv(dri_nutrients, file = "data/clean-data/02_DietaryReferenceIntakes.csv",
  row.names = FALSE)




