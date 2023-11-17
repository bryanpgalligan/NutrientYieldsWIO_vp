## This script estimates WIO nutrient losses due to overfishing only

# Load catch compositions
wio <- readRDS("data/clean-data/04_NationalCatchCompositions.rds")

# Load species data
species <- read_csv("data/clean-data/02_SpeciesData.csv")

# Load DRIs
dri <- read_csv("data/clean-data/02_DietaryReferenceIntakes.csv")

# Load model estimates
reef_cells <- read_csv("data/raw-data/wcs-biomass-model/wcs-biomass-model.csv")

# Make countries lower case
reef_cells$country <- str_to_lower(reef_cells$country)

# Function to find SE (standard error) of the mean of a vector
se.mean <- function(x) sd(x)/sqrt(length(x))




##### Calculate Nutrient Losses #####

## Make a data frame summarizing each country

# Create initial data frame
wio_nutrients <- data.frame(names(wio))
colnames(wio_nutrients) <- "eez"

# People and biomass columns
wio_nutrients$lost.yield_total.tons.day <- NA
wio_nutrients$lost.yield_mean.kg.km2.day <- NA
wio_nutrients$lost.yield_se.mean <- NA
wio_nutrients$population <- c( #source is world bank
  836770, #comoros
  54027490, #kenya
  29611710, #madagascar
  523790, #maldives
  1262520, #mauritius
  310022, #mayotte (source is INSEE)
  32969520, #mozambique
  873102, #reunion (source is INSEE)
  100060, #seychelles
  65497750 #tanzania
  )
wio_nutrients$children <- c( 
  #estimated number of children 1-3
  #source is United Nations, Department of Economic and Social Affairs, Population Division (2022). World Population Prospects: The 2022 Revision, custom data acquired via website.
  68660, #comoros
  4207317, #kenya
  2562401, #madagascar
  21821, #maldives
  39560, #mauritius
  32023, #mayotte
  3274842, #mozambique
  42090, #reunion
  4915, #seychelles
  5624229 #tanzania
)

# Nutrient density columns
wio_nutrients$calcium_pdv.100g <- NA
wio_nutrients$iron_pdv.100g <- NA
wio_nutrients$omega3_pdv.100g <- NA
wio_nutrients$selenium_pdv.100g <- NA
wio_nutrients$vitamina_pdv.100g <- NA
wio_nutrients$zinc_pdv.100g <- NA

# Nutrient yield columns
wio_nutrients$lost.calcium_dri.day <- NA
wio_nutrients$lost.calcium_dri.child.day <- NA
wio_nutrients$lost.calcium_mean.dri.km2.day <- NA
wio_nutrients$lost.calcium_se.dri <- NA

wio_nutrients$lost.iron_dri.day <- NA
wio_nutrients$lost.iron_dri.child.day <- NA
wio_nutrients$lost.iron_mean.dri.km2.day <- NA
wio_nutrients$lost.iron_se.dri <- NA

wio_nutrients$lost.omega3_dri.day <- NA
wio_nutrients$lost.omega3_dri.child.day <- NA
wio_nutrients$lost.omega3_mean.dri.km2.day <- NA
wio_nutrients$lost.omega3_se.dri <- NA

wio_nutrients$lost.selenium_dri.day <- NA
wio_nutrients$lost.selenium_dri.child.day <- NA
wio_nutrients$lost.selenium_mean.dri.km2.day <- NA
wio_nutrients$lost.selenium_se.dri <- NA

wio_nutrients$lost.vitamina_dri.day <- NA
wio_nutrients$lost.vitamina_dri.child.day <- NA
wio_nutrients$lost.vitamina_mean.dri.km2.day <- NA
wio_nutrients$lost.vitamina_se.dri <- NA

wio_nutrients$lost.zinc_dri.day <- NA
wio_nutrients$lost.zinc_dri.child.day <- NA
wio_nutrients$lost.zinc_mean.dri.km2.day <- NA
wio_nutrients$lost.zinc_se.dri <- NA

# Populate data frame with nutrient composition info
for(i in 1:nrow(wio_nutrients)) {
  
  # Get the data frame in question
  eez <- data.frame(wio[[i]])
  
  # Add new columns
  eez$calcium_pdv <- NA
  eez$iron_pdv <- NA
  eez$omega3_pdv <- NA
  eez$selenium_pdv <- NA
  eez$vitamina_pdv <- NA
  eez$zinc_pdv <- NA
  
  # Populate columns with species' nutrient concentrations
  for(j in 1:nrow(eez)){
    
    # If this taxon is a species
    if(eez$rank[j] == "Species"){
      
      # Does this species have nutrient predictions?
      if(is.na(species$calcium_mg.100g[species$species == eez$valid_name[j]]) == FALSE){
        
        # Populate concentration columns
        eez$calcium_pdv[j] <- species$calcium_mg.100g[species$species == eez$valid_name[j]] / dri[dri$nutrients == "calcium_mg", "dri"]
        eez$iron_pdv[j] <- species$iron_mg.100g[species$species == eez$valid_name[j]] / dri[dri$nutrients == "iron_mg", "dri"]
        eez$omega3_pdv[j] <- species$omega3_g.100g[species$species == eez$valid_name[j]] / dri[dri$nutrients == "omega3_g", "dri"]
        eez$selenium_pdv[j] <- species$selenium_ug.100g[species$species == eez$valid_name[j]] / dri[dri$nutrients == "selenium_ug", "dri"]
        eez$vitamina_pdv[j] <- species$vitamin.a_ug.100g[species$species == eez$valid_name[j]] / dri[dri$nutrients == "vitamina_ug", "dri"]
        eez$zinc_pdv[j] <- species$zinc_mg.100g[species$species == eez$valid_name[j]] / dri[dri$nutrients == "zinc_mg", "dri"]
        
      } else{ # If there is no prediction for this species
        
        # Change rank to genus
        eez$rank[j] <- "Genus"
        
      }
      
    }
    
    # If this taxon needs a genus level estimate
    if(eez$rank[j] == "Genus"){
      
      # Subset this genus from species
      species2 <- species[species$genus == eez$genus[j],]
      
      # If the genus has nutrient estimates
      if(is.nan(mean(species2$calcium_mg.100g, na.rm = TRUE)) == FALSE){
        
        # Populate concentration columns
        eez$calcium_pdv[j] <- mean(species2$calcium_mg.100g, na.rm = TRUE) / dri[dri$nutrients == "calcium_mg", "dri"]
        eez$iron_pdv[j] <- mean(species2$iron_mg.100g, na.rm = TRUE) / dri[dri$nutrients == "iron_mg", "dri"]
        eez$omega3_pdv[j] <- mean(species2$omega3_g.100g, na.rm = TRUE) / dri[dri$nutrients == "omega3_g", "dri"]
        eez$selenium_pdv[j] <- mean(species2$selenium_ug.100g, na.rm = TRUE) / dri[dri$nutrients == "selenium_ug", "dri"]
        eez$vitamina_pdv[j] <- mean(species2$vitamin.a_ug.100g, na.rm = TRUE) / dri[dri$nutrients == "vitamina_ug", "dri"]
        eez$zinc_pdv[j] <- mean(species2$zinc_mg.100g, na.rm = TRUE) / dri[dri$nutrients == "zinc_mg", "dri"]
        
      } else{
        
        # If not, assign rank as family
        eez$rank[j] <- "Family"
        
      }
      
    }
    
    # If taxon needs a family level estimate
    if(eez$rank[j] == "Family"){
      
      # Subset species data for this family
      species2 <- species[species$family == eez$family[j], ]
      
      # Populate concentration columns
      eez$calcium_pdv[j] <- mean(species2$calcium_mg.100g, na.rm = TRUE) / dri[dri$nutrients == "calcium_mg", "dri"]
      eez$iron_pdv[j] <- mean(species2$iron_mg.100g, na.rm = TRUE) / dri[dri$nutrients == "iron_mg", "dri"]
      eez$omega3_pdv[j] <- mean(species2$omega3_g.100g, na.rm = TRUE) / dri[dri$nutrients == "omega3_g", "dri"]
      eez$selenium_pdv[j] <- mean(species2$selenium_ug.100g, na.rm = TRUE) / dri[dri$nutrients == "selenium_ug", "dri"]
      eez$vitamina_pdv[j] <- mean(species2$vitamin.a_ug.100g, na.rm = TRUE) / dri[dri$nutrients == "vitamina_ug", "dri"]
      eez$zinc_pdv[j] <- mean(species2$zinc_mg.100g, na.rm = TRUE) / dri[dri$nutrients == "zinc_mg", "dri"]
        
    }
    
    # If taxon needs an order level estimate
    if(eez$rank[j] == "Order"){
      
      # Subset species data for this family
      species2 <- species[species$order == eez$order[j], ]
      
      # Populate concentration columns
      eez$calcium_pdv[j] <- mean(species2$calcium_mg.100g, na.rm = TRUE) / dri[dri$nutrients == "calcium_mg", "dri"]
      eez$iron_pdv[j] <- mean(species2$iron_mg.100g, na.rm = TRUE) / dri[dri$nutrients == "iron_mg", "dri"]
      eez$omega3_pdv[j] <- mean(species2$omega3_g.100g, na.rm = TRUE) / dri[dri$nutrients == "omega3_g", "dri"]
      eez$selenium_pdv[j] <- mean(species2$selenium_ug.100g, na.rm = TRUE) / dri[dri$nutrients == "selenium_ug", "dri"]
      eez$vitamina_pdv[j] <- mean(species2$vitamin.a_ug.100g, na.rm = TRUE) / dri[dri$nutrients == "vitamina_ug", "dri"]
      eez$zinc_pdv[j] <- mean(species2$zinc_mg.100g, na.rm = TRUE) / dri[dri$nutrients == "zinc_mg", "dri"]
        
      
    }
    
    # If taxon needs a class level estimate
    if(eez$rank[j] == "Class"){
      
      # Subset species data for this family
      species2 <- species[species$class == eez$class[j], ]
      
      # Populate concentration columns
      eez$calcium_pdv[j] <- mean(species2$calcium_mg.100g, na.rm = TRUE) / dri[dri$nutrients == "calcium_mg", "dri"]
      eez$iron_pdv[j] <- mean(species2$iron_mg.100g, na.rm = TRUE) / dri[dri$nutrients == "iron_mg", "dri"]
      eez$omega3_pdv[j] <- mean(species2$omega3_g.100g, na.rm = TRUE) / dri[dri$nutrients == "omega3_g", "dri"]
      eez$selenium_pdv[j] <- mean(species2$selenium_ug.100g, na.rm = TRUE) / dri[dri$nutrients == "selenium_ug", "dri"]
      eez$vitamina_pdv[j] <- mean(species2$vitamin.a_ug.100g, na.rm = TRUE) / dri[dri$nutrients == "vitamina_ug", "dri"]
      eez$zinc_pdv[j] <- mean(species2$zinc_mg.100g, na.rm = TRUE) / dri[dri$nutrients == "zinc_mg", "dri"]
        
    }
    
  }
  
  # Make nutrient columns not lists
  eez$calcium_pdv <- unlist(eez$calcium_pdv)
  eez$iron_pdv <- unlist(eez$iron_pdv)
  eez$omega3_pdv <- unlist(eez$omega3_pdv)
  eez$selenium_pdv <- unlist(eez$selenium_pdv)
  eez$vitamina_pdv <- unlist(eez$vitamina_pdv)
  eez$zinc_pdv <- unlist(eez$zinc_pdv)
  
  # Find weighted means and populate wio_nutrients
  wio_nutrients$calcium_pdv.100g[i] <- weighted.mean(eez$calcium_pdv, eez$ratio, na.rm = TRUE)
  wio_nutrients$iron_pdv.100g[i] <- weighted.mean(eez$iron_pdv, eez$ratio, na.rm = TRUE)
  wio_nutrients$omega3_pdv.100g[i] <- weighted.mean(eez$omega3_pdv, eez$ratio, na.rm = TRUE)
  wio_nutrients$selenium_pdv.100g[i] <- weighted.mean(eez$selenium_pdv, eez$ratio, na.rm = TRUE)
  wio_nutrients$vitamina_pdv.100g[i] <- weighted.mean(eez$vitamina_pdv, eez$ratio, na.rm = TRUE)
  wio_nutrients$zinc_pdv.100g[i] <- weighted.mean(eez$zinc_pdv, eez$ratio, na.rm = TRUE)
  
}


# Populate lost yields based on biomass estimates, only including losses due to overfishing
for (i in 1:nrow(wio_nutrients)) {
  
  # Subset model estimates for this eez
  eez <- filter(reef_cells, country == wio_nutrients$eez[i])
  
  # Set all underfished and sustainable rows to 0
  eez$lost.yield[eez$b.bmsy > 1] <- 0
  
  # Subset only one row for each cell
  eez <- filter(eez, vertex.index == 1)
  
  # Find the total lost yield for each cell
  eez$total.lost.yield <- eez$lost.yield * eez$area_km2
  
  # Find total lost yield for eez
  wio_nutrients$lost.yield_total.tons.day[i] <- sum(eez$total.lost.yield) / 365
  
  # Find mean lost yield per km2 per day
  a <- mean(eez$lost.yield)
  a <- a * 1000 #convert to kg
  a <- a / 365 #convert to daily
  wio_nutrients$lost.yield_mean.kg.km2.day[i] <- a
  
  # Find standard error of mean for lost yield
  b <- se.mean(eez$lost.yield)
  b <- b * 1000 #convert to kg
  b <- b / 365 #convert to daily
  wio_nutrients$lost.yield_se.mean[i] <- b
  
}

# Calculate lost nutrients by eez
wio_nutrients$lost.calcium_dri.day <- (wio_nutrients$lost.yield_total.tons.day * 10000) * wio_nutrients$calcium_pdv.100g
wio_nutrients$lost.calcium_dri.child.day <- wio_nutrients$lost.calcium_dri.day / wio_nutrients$children
wio_nutrients$lost.calcium_mean.dri.km2.day <- (wio_nutrients$lost.yield_mean.kg.km2.day * 10) * wio_nutrients$calcium_pdv.100g
wio_nutrients$lost.calcium_se.dri <- (wio_nutrients$lost.yield_se.mean * 10) * wio_nutrients$calcium_pdv.100g

wio_nutrients$lost.iron_dri.day <- (wio_nutrients$lost.yield_total.tons.day * 10000) * wio_nutrients$iron_pdv.100g
wio_nutrients$lost.iron_dri.child.day <- wio_nutrients$lost.iron_dri.day / wio_nutrients$children
wio_nutrients$lost.iron_mean.dri.km2.day <- (wio_nutrients$lost.yield_mean.kg.km2.day * 10) * wio_nutrients$iron_pdv.100g
wio_nutrients$lost.iron_se.dri <- (wio_nutrients$lost.yield_se.mean * 10) * wio_nutrients$iron_pdv.100g

wio_nutrients$lost.omega3_dri.day <- (wio_nutrients$lost.yield_total.tons.day * 10000) * wio_nutrients$omega3_pdv.100g
wio_nutrients$lost.omega3_dri.child.day <- wio_nutrients$lost.omega3_dri.day / wio_nutrients$children
wio_nutrients$lost.omega3_mean.dri.km2.day <- (wio_nutrients$lost.yield_mean.kg.km2.day * 10) * wio_nutrients$omega3_pdv.100g
wio_nutrients$lost.omega3_se.dri <- (wio_nutrients$lost.yield_se.mean * 10) * wio_nutrients$omega3_pdv.100g

wio_nutrients$lost.selenium_dri.day <- (wio_nutrients$lost.yield_total.tons.day * 10000) * wio_nutrients$selenium_pdv.100g
wio_nutrients$lost.selenium_dri.child.day <- wio_nutrients$lost.selenium_dri.day / wio_nutrients$children
wio_nutrients$lost.selenium_mean.dri.km2.day <- (wio_nutrients$lost.yield_mean.kg.km2.day * 10) * wio_nutrients$selenium_pdv.100g
wio_nutrients$lost.selenium_se.dri <- (wio_nutrients$lost.yield_se.mean * 10) * wio_nutrients$selenium_pdv.100g

wio_nutrients$lost.vitamina_dri.day <- (wio_nutrients$lost.yield_total.tons.day * 10000) * wio_nutrients$vitamina_pdv.100g
wio_nutrients$lost.vitamina_dri.child.day <- wio_nutrients$lost.vitamina_dri.day / wio_nutrients$children
wio_nutrients$lost.vitamina_mean.dri.km2.day <- (wio_nutrients$lost.yield_mean.kg.km2.day * 10) * wio_nutrients$vitamina_pdv.100g
wio_nutrients$lost.vitamina_se.dri <- (wio_nutrients$lost.yield_se.mean * 10) * wio_nutrients$vitamina_pdv.100g

wio_nutrients$lost.zinc_dri.day <- (wio_nutrients$lost.yield_total.tons.day * 10000) * wio_nutrients$zinc_pdv.100g
wio_nutrients$lost.zinc_dri.child.day <- wio_nutrients$lost.zinc_dri.day / wio_nutrients$children
wio_nutrients$lost.zinc_mean.dri.km2.day <- (wio_nutrients$lost.yield_mean.kg.km2.day * 10) * wio_nutrients$zinc_pdv.100g
wio_nutrients$lost.zinc_se.dri <- (wio_nutrients$lost.yield_se.mean * 10) * wio_nutrients$zinc_pdv.100g


# Add a total nutrients column - gives you nutrient density out of 1
wio_nutrients$nutrient.density_coef <- (wio_nutrients$calcium_pdv.100g + 
  wio_nutrients$iron_pdv.100g +
  wio_nutrients$omega3_pdv.100g + 
  #wio_nutrients$selenium_pdv.100g +
  1 + # replace selenium with 1 because all values are >1
  wio_nutrients$vitamina_pdv.100g +
  wio_nutrients$zinc_pdv.100g) / 6

# Add total lost nutrients and SE columns
wio_nutrients$lost.dri_km2.day <- (wio_nutrients$lost.yield_mean.kg.km2.day * 10) * wio_nutrients$nutrient.density_coef
wio_nutrients$lost.dri_se <- (wio_nutrients$lost.yield_se.mean * 10) * wio_nutrients$nutrient.density_coef

# Save WIO nutrients
write.csv(wio_nutrients, file = "data/clean-data/05_NationalNutrientLosses.csv", row.names = FALSE)


# Remove unused countries
reef_cells <- filter(reef_cells, country %in% wio_nutrients$eez)

# Save a second version of reef_points
reef_points <- filter(reef_cells, vertex.index == 1)

# Remove all vertex indices that = 0
reef_cells <- filter(reef_cells, vertex.index != 0)

# Make all underfished sites = 0
reef_cells$lost.yield[reef_cells$b.bmsy > 1] <- 0

# Make lost yields kg per km2 per day
reef_cells$lost.yield <- (reef_cells$lost.yield * 1000) / 365


# Add lost nutrients to model estimates
reef_cells$lost.nutrients <- NA

# Loop to fill in empty column
for(i in 1:nrow(reef_cells)){
  
  # Extract nutrient density for this country
  a <- which(wio_nutrients$eez == reef_cells$country[i])
  a <- wio_nutrients$nutrient.density_coef[a]
  
  # Convert kg/day into 100g portions/day
  b <- reef_cells$lost.yield[i] * 10
  
  # Calculate lost nutrients
  c <- b*a
  
  # Save lost nutrient yield
  reef_cells$lost.nutrients[i] <- c
  
}




##### Nutrient Densities #####

# Data frame of nutrient compositions
temp <- wio_nutrients

# Limit selenium to 100%
temp$selenium_pdv.100g[temp$selenium_pdv.100g > 1] <- 1

# Add a column for total nutrients
temp$total_pdv.100g <- NA
for(i in 1:nrow(temp)){
  
  temp$total_pdv.100g[i] <- sum(temp$calcium_pdv.100g[i], temp$iron_pdv.100g[i],
    temp$omega3_pdv.100g[i], temp$selenium_pdv.100g[i], temp$vitamina_pdv.100g[i],
    temp$zinc_pdv.100g[i])
  
}

# Reformat data
temp <- pivot_longer(temp, cols = calcium_pdv.100g:zinc_pdv.100g, names_to = "nutrients")

# Reorder nutrients
temp$nutrients <- fct_reorder(temp$nutrients, temp$value)

# Make EEZs capital
temp$eez <- str_to_title(temp$eez)

# Reorder EEZs
temp$eez <- fct_reorder(temp$eez, temp$total_pdv.100g)

# Select important data only
temp <- select(temp, eez, lost.yield_mean.kg.km2.day, nutrients, value)

# Plot data for fig 1
ggplot(temp, aes(x = eez, y = value, fill = nutrients, label = (round(value * 100)))) +
  geom_bar(position = "stack", stat = "identity") +
  geom_text(size = 2, position = position_stack(vjust = 0.5), color = "white") +
  scale_y_continuous(labels = percent) +
  ylab(expression("Nutrients in 100g (% Recommended Intake)")) +
  xlab("") +
  coord_flip() +
  scale_fill_jco(
    breaks = c("selenium_pdv.100g", "zinc_pdv.100g", "omega3_pdv.100g", "iron_pdv.100g",
      "vitamina_pdv.100g", "calcium_pdv.100g"),
    labels = c("Selenium", "Omega-3", "Zinc", "Iron", "Vitamin A", "Calcium")) +
  guides(fill = guide_legend(title = "", nrow = 1, label.position = "bottom"
    #keywidth = unit(1, "mm"), keyheight = unit(1, "mm")
    )) +
  theme_pubr() +
  theme(legend.position = "bottom",
    text = element_text(size = 6.5),
    legend.text = element_text(size = 6.5),
    axis.title.x = element_text(vjust = -5),
    plot.margin = unit(c(1, 5, 1, 1), "mm"),
    legend.box.spacing = unit(3, "mm")
    )

# Save fig 1
ggsave("figures/05_NationalNutrientDensities.pdf",
  width = 114, height = 85, units = "mm",
  dpi = 300)




##### Heat Maps of Losses #####

## Reunion

# reunion reef cells
reunion_reefs <- filter(reef_cells, country == "reunion")

# Get map data
reunion_map <- get_stamenmap(maptype = "terrain-background",
  bbox = c(left = min(reunion_reefs$xcoord) - 0.1,
    bottom = min(reunion_reefs$ycoord) - 0.1,
    right = max(reunion_reefs$xcoord) + 0.1,
    top = max(reunion_reefs$ycoord) + 0.1))

# Plot map
reunion <- ggmap(reunion_map) +
  geom_polygon(data = reunion_reefs,
    aes(x = xcoord, y = ycoord,
      fill = lost.nutrients, group = id,
      linetype = NA)) +
  ggtitle("Reunion") +
  xlab("") +
  ylab("") +
  theme_pubr() +
  scale_fill_gradient(low = "blue", high = "red", limits = c(0, 31),
    name = expression(paste(
      "Lost Nutrient \nYield ",
      bgroup("(", frac("DRI", km^2~day), ")"),
      ))) +
  #guides(fill = guide_colorbar(title.position = "top")) +
  theme(text = element_text(size = 6.5), title = element_text(face = "bold"),
    legend.position = "right", legend.justification = "center",
    legend.title = element_text(vjust = 0)
    )


## NE Mozambique

# mozambique reef cells
mozambique_reefs <- filter(reef_cells, country == "mozambique")

# Get map data
mozambique_map <- get_stamenmap(maptype = "terrain-background",
  bbox = c(left = 39,
    bottom = -16,
    right = 42,
    top = max(mozambique_reefs$ycoord) + 0.1))

# Plot map
mozambique <- ggmap(mozambique_map) +
  geom_polygon(data = mozambique_reefs,
    aes(x = xcoord, y = ycoord,
      fill = lost.nutrients, group = id,
      linetype = NA)) +
  ggtitle("Northeast Mozambique") +
  xlab("") +
  ylab("") +
  theme_pubr() +
  scale_fill_gradient(low = "blue", high = "red", limits = c(0, 31),
    name = expression(paste("Daily Lost Nutrients ", bgroup("(", frac("RDI", km^2), ")")))) +
  theme(text = element_text(size = 6.5), title = element_text(face = "bold"),
    legend.position = "none")


## Tanzania

# tanzania reef cells
tanzania_reefs <- filter(reef_cells, country == "tanzania")

# Get map data
tanzania_map <- get_stamenmap(maptype = "terrain-background",
  bbox = c(left = 37.5,
    bottom = min(tanzania_reefs$ycoord) - 0.1,
    right = max(tanzania_reefs$xcoord) + 0.1,
    top = max(tanzania_reefs$ycoord) + 0.1))

# Plot map
tanzania <- ggmap(tanzania_map) +
  geom_polygon(data = tanzania_reefs,
    aes(x = xcoord, y = ycoord,
      fill = lost.nutrients, group = id,
      linetype = NA)) +
  ggtitle("Tanzania") +
  xlab("") +
  ylab("") +
  scale_x_continuous(breaks = c(38, 39, 40, 41)) +
  theme_pubr() +
  scale_fill_gradient(low = "blue", high = "red", limits = c(0, 31),
    name = expression(paste("Daily Lost Nutrients ", bgroup("(", frac("RDI", km^2), ")")))) +
  theme(text = element_text(size = 6.5), title = element_text(face = "bold"),
    legend.position = "none")


## Kenya

# kenya reef cells
kenya_reefs <- filter(reef_cells, country == "kenya")

# Get map data
kenya_map <- get_stamenmap(maptype = "terrain-background",
  bbox = c(left = min(kenya_reefs$xcoord) - 0.1,
    bottom = min(kenya_reefs$ycoord) - 0.1,
    right = max(kenya_reefs$xcoord) + 0.1,
    top = max(kenya_reefs$ycoord) + 0.1))

# Plot map
kenya <- ggmap(kenya_map) +
  geom_polygon(data = kenya_reefs,
    aes(x = xcoord, y = ycoord,
      fill = lost.nutrients, group = id,
      linetype = NA)) +
  ggtitle("Kenya") +
  xlab("") +
  ylab("") +
  theme_pubr() +
  scale_fill_gradient(low = "blue", high = "red", limits = c(0, 31),
    name = expression(paste("Daily Lost Nutrients ", bgroup("(", frac("RDI", km^2), ")")))) +
  theme(text = element_text(size = 6.5), title = element_text(face = "bold"),
    legend.position = "none")


## Mauritius

# mauritius reef cells
mauritius_reefs <- filter(reef_cells, country == "mauritius")

# Get map data
mauritius_map <- get_stamenmap(maptype = "terrain-background",
  bbox = c(left = 57,
    bottom = min(mauritius_reefs$ycoord) - 0.1,
    right = 58,
    top = -19.7))

# Plot map
mauritius <- ggmap(mauritius_map) +
  geom_polygon(data = mauritius_reefs,
    aes(x = xcoord, y = ycoord,
      fill = lost.nutrients, group = id,
      linetype = NA)) +
  ggtitle("Mauritius") +
  xlab("") +
  ylab("") +
  #scale_x_continuous(breaks = c(57, 57.5, 58)) +
  theme_pubr() +
  scale_fill_gradient(low = "blue", high = "red", limits = c(0, 31),
    name = expression(paste("Daily Lost Nutrients ", bgroup("(", frac("RDI", km^2), ")")))) +
  theme(text = element_text(size = 6.5), title = element_text(face = "bold"),
    legend.position = "none",
    plot.margin = unit(c(0, 4, 0, 0), "mm"))


## SW Madagascar

# madagascar reef cells
madagascar_reefs <- filter(reef_cells, country == "madagascar")

# Get map data
madagascar_map <- get_stamenmap(maptype = "terrain-background",
  bbox = c(left = min(madagascar_reefs$xcoord) - 0.1,
    bottom = min(madagascar_reefs$ycoord) - 0.1,
    right = 45,
    top = -23))

# Plot map
madagascar <- ggmap(madagascar_map) +
  geom_polygon(data = madagascar_reefs,
    aes(x = xcoord, y = ycoord,
      fill = lost.nutrients, group = id,
      linetype = NA)) +
  ggtitle("Southwest Madagascar") +
  xlab("") +
  ylab("") +
  theme_pubr() +
  scale_fill_gradient(low = "blue", high = "red", limits = c(0, 31),
    name = expression(paste("Daily Lost Nutrients ", bgroup("(", frac("RDI", km^2), ")")))) +
  theme(text = element_text(size = 6.5), title = element_text(face = "bold"),
    legend.position = "none")


# ## Legend
# 
# # Save the reunion legend as a separate ggplot object
# legend_extract <- function(ggplot_obj) {
#   
#   # Convert the ggplot object to a grob
#   ggplot_grob <- ggplotGrob(ggplot_obj)
#   
#   # Find the position of the legend in the grob
#   legend_pos <- which(ggplot_grob$layout$name == "guide-box")
#   
#   # Extract the legend from the grob
#   legend_grob <- ggplot_grob$grobs[[legend_pos]]
#   
#   # Create a new ggplot object with just the legend
#   legend_plot <- ggplot() +
#     annotation_custom(legend_grob)
#   
#   return(legend_plot)
# }
# 
# # Call the function to extract the legend as a separate ggplot object
# leg <- legend_extract(reunion) +
#   theme_pubclean() +
#   theme(text = element_text(size = 6.5),
#     plot.margin = unit(c(0, 0, 0, 20), "mm"))


## Combine plots
ggarrange(reunion, mozambique, kenya,
  tanzania, mauritius, madagascar,
  ncol = 3, nrow = 2,
  common.legend = TRUE, legend.grob = get_legend(reunion), legend = "right")

# Save figure
ggsave("figures/05_ReefMaps.pdf",
  width = 174, units = "mm",
  dpi = 300)


## Calculate average losses for overfished sites
overfished <- filter(reef_cells, status == "Overfished")
mean(overfished$lost.yield, na.rm = TRUE)
se.mean(overfished$lost.yield)
table(reef_cells$status)




##### Supplementary Heat Maps #####

## Reunion

# Get map data
reunion_map <- get_stamenmap(maptype = "terrain",
  bbox = c(left = min(reunion_reefs$xcoord) - 0.1,
    bottom = min(reunion_reefs$ycoord) - 0.1,
    right = 55.9,
    top = -20.8),
  force = TRUE)

# # Get labels
# reunion_labels <- get_stamenmap(maptype = "terrain-labels",
#   bbox = c(left = min(reunion_reefs$xcoord) - 0.1,
#     bottom = min(reunion_reefs$ycoord) - 0.1,
#     right = 55.9,
#     top = -20.8),
#   zoom = 12,
#   force = TRUE)

# Plot map
ggmap(reunion_map) +
  geom_polygon(data = reunion_reefs,
    aes(x = xcoord, y = ycoord,
      fill = lost.nutrients, group = id,
      linetype = NA)) +
  # inset_ggmap(reunion_labels) +
  xlab("") +
  ylab("") +
  theme_pubr() +
  scale_fill_gradient(low = "blue", high = "red", limits = c(0, 31),
    name = expression(paste(
      "Lost Nutrient Yield ",
      bgroup("(", frac("DRI", km^2~day), ")"),
      ))) +
  #guides(fill = guide_colorbar(title.position = "top")) +
  theme(text = element_text(size = 7.5), title = element_text(face = "bold"),
    legend.position = "bottom", legend.justification = "center",
    legend.title = element_text(vjust = 0),
    plot.margin = unit(c(0, 5, 0, 0), units = "mm")
    )

# Save image
ggsave("figures/05_ReunionMap.pdf",
  width = 174, height = 200, units = "mm",
  dpi = 300)


## Mozambique

# Get map data
mozambique_map <- get_stamenmap(maptype = "terrain-background",
  bbox = c(left = min(mozambique_reefs$xcoord) - 0.1,
    bottom = min(mozambique_reefs$ycoord) - 0.1,
    right = 42,
    top = max(mozambique_reefs$ycoord) + 0.1),
  zoom = 7, crop = FALSE)

# Get labels
mozambique_labels <- get_stamenmap(maptype = "terrain-labels",
  bbox = c(left = min(mozambique_reefs$xcoord) - 0.1,
    bottom = min(mozambique_reefs$ycoord) - 0.1,
    right = 42,
    top = max(mozambique_reefs$ycoord) + 0.1),
  zoom = 7, crop = FALSE)

# Plot map
ggmap(mozambique_map) +
  geom_polygon(data = mozambique_reefs,
    aes(x = xcoord, y = ycoord,
      fill = lost.nutrients, group = id,
      linetype = NA)) +
  inset_ggmap(mozambique_labels) +
  xlab("") +
  ylab("") +
  theme_pubr() +
  scale_fill_gradient(low = "blue", high = "red", limits = c(0, 31),
    name = expression(paste(
      "Lost Nutrient Yield ",
      bgroup("(", frac("DRI", km^2~day), ")"),
      ))) +
  theme(text = element_text(size = 7.5), title = element_text(face = "bold"),
    legend.position = "bottom", legend.justification = "center",
    legend.title = element_text(vjust = 0)
    )

# Save image
ggsave("figures/05_MozambiqueMap.pdf",
  width = 174, height = 200, units = "mm",
  dpi = 300)


## Tanzania

# Get map data
tanzania_map <- get_stamenmap(maptype = "terrain-background",
  bbox = c(left = 37.5,
    bottom = min(tanzania_reefs$ycoord) - 0.1,
    right = max(tanzania_reefs$xcoord) + 0.1,
    top = max(tanzania_reefs$ycoord) + 0.1),
  zoom = 7)

# Get labels
tanzania_labels <- get_stamenmap(maptype = "terrain-labels",
  bbox = c(left = 37.5,
    bottom = min(tanzania_reefs$ycoord) - 0.1,
    right = max(tanzania_reefs$xcoord) + 0.1,
    top = max(tanzania_reefs$ycoord) + 0.1),
  zoom = 8)

# Plot map
ggmap(tanzania_map) +
  geom_polygon(data = tanzania_reefs,
    aes(x = xcoord, y = ycoord,
      fill = lost.nutrients, group = id,
      linetype = NA)) +
  inset_ggmap(tanzania_labels) +
  xlab("") +
  ylab("") +
  theme_pubr() +
  scale_fill_gradient(low = "blue", high = "red", limits = c(0, 31),
    name = expression(paste(
      "Lost Nutrient Yield ",
      bgroup("(", frac("DRI", km^2~day), ")"),
      ))) +
  theme(text = element_text(size = 7.5), title = element_text(face = "bold"),
    legend.position = "bottom", legend.justification = "center",
    legend.title = element_text(vjust = 0)
    )

# Save image
ggsave("figures/05_TanzaniaMap.pdf",
  width = 174, height = 200, units = "mm",
  dpi = 300)


## Kenya

# Get map data
kenya_map <- get_stamenmap(maptype = "terrain-background",
  bbox = c(left = min(kenya_reefs$xcoord) - 0.1,
    bottom = min(kenya_reefs$ycoord) - 0.1,
    right = max(kenya_reefs$xcoord) + 0.1,
    top = max(kenya_reefs$ycoord) + 0.1),
  zoom = 8)

# Get labels
kenya_labels <- get_stamenmap(maptype = "terrain-labels",
  bbox = c(left = min(kenya_reefs$xcoord) - 0.1,
    bottom = min(kenya_reefs$ycoord) - 0.1,
    right = max(kenya_reefs$xcoord) + 0.1,
    top = max(kenya_reefs$ycoord) + 0.1),
  zoom = 8)

# Plot map
ggmap(kenya_map) +
  geom_polygon(data = kenya_reefs,
    aes(x = xcoord, y = ycoord,
      fill = lost.nutrients, group = id,
      linetype = NA)) +
  inset_ggmap(kenya_labels) +
  xlab("") +
  ylab("") +
  theme_pubr() +
  scale_fill_gradient(low = "blue", high = "red", limits = c(0, 31),
    name = expression(paste(
      "Lost Nutrient Yield ",
      bgroup("(", frac("DRI", km^2~day), ")"),
      ))) +
  theme(text = element_text(size = 7.5), title = element_text(face = "bold"),
    legend.position = "bottom", legend.justification = "center",
    legend.title = element_text(vjust = 0)
    )

# Save image
ggsave("figures/05_KenyaMap.pdf",
  width = 174, height = 200, units = "mm",
  dpi = 300)


## Mauritius

# Get map data
mauritius_map <- get_stamenmap(maptype = "terrain-background",
  bbox = c(left = min(mauritius_reefs$xcoord) - 0.1,
    bottom = min(mauritius_reefs$ycoord) - 0.1,
    right = max(mauritius_reefs$xcoord) + 0.1,
    top = max(mauritius_reefs$ycoord) + 0.1),
  zoom = 8)

# # Get labels
# mauritius_labels <- get_stamenmap(maptype = "terrain-labels",
#   bbox = c(left = min(mauritius_reefs$xcoord) - 0.1,
#     bottom = min(mauritius_reefs$ycoord) - 0.1,
#     right = max(mauritius_reefs$xcoord) + 0.1,
#     top = max(mauritius_reefs$ycoord) + 0.1),
#   zoom = 8)

# Plot map
ggmap(mauritius_map) +
  geom_polygon(data = mauritius_reefs,
    aes(x = xcoord, y = ycoord,
      fill = lost.nutrients, group = id,
      linetype = NA)) +
  #inset_ggmap(mauritius_labels) +
  xlab("") +
  ylab("") +
  theme_pubr() +
  scale_fill_gradient(low = "blue", high = "red", limits = c(0, 31),
    name = expression(paste(
      "Lost Nutrient Yield ",
      bgroup("(", frac("DRI", km^2~day), ")"),
      ))) +
  theme(text = element_text(size = 7.5), title = element_text(face = "bold"),
    legend.position = "bottom", legend.justification = "center",
    legend.title = element_text(vjust = 0)
    )

# Save image
ggsave("figures/05_MauritiusMap.pdf",
  width = 174, height = 200, units = "mm",
  dpi = 300)


## Madagascar

# Get map data
madagascar_map <- get_stamenmap(maptype = "terrain-background",
  bbox = c(left = min(madagascar_reefs$xcoord) - 0.1,
    bottom = min(madagascar_reefs$ycoord) - 0.3,
    right = max(madagascar_reefs$xcoord) + 0.1,
    top = max(madagascar_reefs$ycoord) + 0.1),
  zoom = 7)

# Get labels
madagascar_labels <- get_stamenmap(maptype = "terrain-labels",
  bbox = c(left = min(madagascar_reefs$xcoord) - 0.1,
    bottom = min(madagascar_reefs$ycoord) - 0.3,
    right = max(madagascar_reefs$xcoord) + 0.1,
    top = max(madagascar_reefs$ycoord) + 0.1),
  zoom = 7)

# Plot map
ggmap(madagascar_map) +
  geom_polygon(data = madagascar_reefs,
    aes(x = xcoord, y = ycoord,
      fill = lost.nutrients, group = id,
      linetype = NA)) +
  inset_ggmap(madagascar_labels) +
  xlab("") +
  ylab("") +
  theme_pubr() +
  scale_fill_gradient(low = "blue", high = "red", limits = c(0, 31),
    name = expression(paste(
      "Lost Nutrient Yield ",
      bgroup("(", frac("DRI", km^2~day), ")"),
      ))) +
  theme(text = element_text(size = 7.5), title = element_text(face = "bold"),
    legend.position = "bottom", legend.justification = "center",
    legend.title = element_text(vjust = 0)
    )

# Save image
ggsave("figures/05_MadagascarMap.pdf",
  width = 174, height = 200, units = "mm",
  dpi = 300)


## Seychelles

# Get reef cells
seychelles_reefs <- filter(reef_cells, country == "seychelles")

# Get map data
seychelles_map <- get_stamenmap(maptype = "terrain-background",
  bbox = c(left = min(seychelles_reefs$xcoord) - 0.4,
    bottom = min(seychelles_reefs$ycoord) - 0.1,
    right = max(seychelles_reefs$xcoord) + 0.1,
    top = max(seychelles_reefs$ycoord) + 0.1),
  zoom = 7)

# Get labels
seychelles_labels <- get_stamenmap(maptype = "terrain-labels",
  bbox = c(left = min(seychelles_reefs$xcoord) - 0.4,
    bottom = min(seychelles_reefs$ycoord) - 0.1,
    right = max(seychelles_reefs$xcoord) + 0.1,
    top = max(seychelles_reefs$ycoord) + 0.1),
  zoom = 7)

# Plot map
ggmap(seychelles_map) +
  geom_polygon(data = seychelles_reefs,
    aes(x = xcoord, y = ycoord,
      fill = lost.nutrients, group = id,
      linetype = NA)) +
  inset_ggmap(seychelles_labels) +
  xlab("") +
  ylab("") +
  theme_pubr() +
  scale_fill_gradient(low = "blue", high = "red", limits = c(0, 31),
    name = expression(paste(
      "Lost Nutrient Yield ",
      bgroup("(", frac("DRI", km^2~day), ")"),
      ))) +
  theme(text = element_text(size = 7.5), title = element_text(face = "bold"),
    legend.position = "bottom", legend.justification = "center",
    legend.title = element_text(vjust = 0)
    )

# Save image
ggsave("figures/05_SeychellesMap.pdf",
  width = 174, height = 200, units = "mm",
  dpi = 300)


## Mayotte

# Get reef cells
mayotte_reefs <- filter(reef_cells, country == "mayotte")

# Get map data
mayotte_map <- get_stamenmap(maptype = "terrain-background",
  bbox = c(left = min(mayotte_reefs$xcoord) - 0.1,
    bottom = min(mayotte_reefs$ycoord) - 0.1,
    right = max(mayotte_reefs$xcoord) + 0.1,
    top = max(mayotte_reefs$ycoord) + 0.1))

# Get labels
mayotte_labels <- get_stamenmap(maptype = "terrain-labels",
  bbox = c(left = min(mayotte_reefs$xcoord) - 0.1,
    bottom = min(mayotte_reefs$ycoord) - 0.1,
    right = max(mayotte_reefs$xcoord) + 0.1,
    top = max(mayotte_reefs$ycoord) + 0.1),
  zoom = 9)

# Plot map
ggmap(mayotte_map) +
  geom_polygon(data = mayotte_reefs,
    aes(x = xcoord, y = ycoord,
      fill = lost.nutrients, group = id,
      linetype = NA)) +
  inset_ggmap(mayotte_labels) +
  xlab("") +
  ylab("") +
  theme_pubr() +
  scale_fill_gradient(low = "blue", high = "red", limits = c(0, 31),
    name = expression(paste(
      "Lost Nutrient Yield ",
      bgroup("(", frac("DRI", km^2~day), ")"),
      ))) +
  theme(text = element_text(size = 7.5), title = element_text(face = "bold"),
    legend.position = "bottom", legend.justification = "center",
    legend.title = element_text(vjust = 0)
    )

# Save image
ggsave("figures/05_MayotteMap.pdf",
  width = 174, height = 200, units = "mm",
  dpi = 300)


## Comoros

# Get reef cells
comoros_reefs <- filter(reef_cells, country == "comoros")

# Get map data
comoros_map <- get_stamenmap(maptype = "terrain-background",
  bbox = c(left = min(comoros_reefs$xcoord) - 0.1,
    bottom = min(comoros_reefs$ycoord) - 0.1,
    right = max(comoros_reefs$xcoord) + 0.1,
    top = max(comoros_reefs$ycoord) + 0.1))

# Get labels
comoros_labels <- get_stamenmap(maptype = "terrain-labels",
  bbox = c(left = min(comoros_reefs$xcoord) - 0.1,
    bottom = min(comoros_reefs$ycoord) - 0.1,
    right = max(comoros_reefs$xcoord) + 0.1,
    top = max(comoros_reefs$ycoord) + 0.1),
  zoom = 9)

# Plot map
ggmap(comoros_map) +
  geom_polygon(data = comoros_reefs,
    aes(x = xcoord, y = ycoord,
      fill = lost.nutrients, group = id,
      linetype = NA)) +
  inset_ggmap(comoros_labels) +
  xlab("") +
  ylab("") +
  theme_pubr() +
  scale_fill_gradient(low = "blue", high = "red", limits = c(0, 31),
    name = expression(paste(
      "Lost Nutrient Yield ",
      bgroup("(", frac("DRI", km^2~day), ")"),
      ))) +
  theme(text = element_text(size = 7.5), title = element_text(face = "bold"),
    legend.position = "bottom", legend.justification = "center",
    legend.title = element_text(vjust = 0)
    )

# Save image
ggsave("figures/05_ComorosMap.pdf",
  width = 174, height = 200, units = "mm",
  dpi = 300)


## Maldives

# Get reef cells
maldives_reefs <- filter(reef_cells, country == "maldives")

# Get map data
maldives_map <- get_stamenmap(maptype = "terrain-background",
  bbox = c(left = min(maldives_reefs$xcoord) - 1,
    bottom = min(maldives_reefs$ycoord) - 0.1,
    right = max(maldives_reefs$xcoord) + 1,
    top = max(maldives_reefs$ycoord) + 0.1),
  zoom = 8)

# Get labels
maldives_labels <- get_stamenmap(maptype = "terrain-labels",
  bbox = c(left = min(maldives_reefs$xcoord) - 1,
    bottom = min(maldives_reefs$ycoord) - 0.1,
    right = max(maldives_reefs$xcoord) + 1,
    top = max(maldives_reefs$ycoord) + 0.1),
  zoom = 7)

# Plot map
ggmap(maldives_map) +
  geom_polygon(data = maldives_reefs,
    aes(x = xcoord, y = ycoord,
      fill = lost.nutrients, group = id,
      linetype = NA)) +
  inset_ggmap(maldives_labels) +
  xlab("") +
  ylab("") +
  theme_pubr() +
  scale_fill_gradient(low = "blue", high = "red", limits = c(0, 31),
    name = expression(paste(
      "Lost Nutrient Yield ",
      bgroup("(", frac("DRI", km^2~day), ")"),
      ))) +
  theme(text = element_text(size = 7.5), title = element_text(face = "bold"),
    legend.position = "bottom", legend.justification = "center",
    legend.title = element_text(vjust = 0)
    )

# Save image
ggsave("figures/05_MaldivesMap.pdf",
  width = 174, height = 200, units = "mm",
  dpi = 300)




##### Plot Nutrient Losses #####

# Make yield kg per km2 per day
reef_points$yield_kg.day <- (reef_points$yield_tons * 1000) / 365

# Make lost yield kg per km2 per day
reef_points$lost.yield_kg.day <- (reef_points$lost.yield * 1000) / 365

# Add lost nutrient columns
reef_points$lost.selenium_dri.day <- NA
reef_points$lost.omega3_dri.day <- NA
reef_points$lost.zinc_dri.day <- NA

# Loop to fill in empty column
for(i in 1:nrow(reef_points)){
  
  # Extract nutrient densities for this country
  a <- which(wio_nutrients$eez == reef_points$country[i])
  se <- 1 #cap selenium density at 100%
  o3 <- wio_nutrients$omega3_pdv.100g[a]
  zn <- wio_nutrients$zinc_pdv.100g[a]
  
  # Convert kg/day into 100g portions/day
  b <- reef_points$lost.yield_kg.day[i] * 10
  
  # Calculate nutrient yields
  reef_points$lost.selenium_dri.day[i] <- b * se
  reef_points$lost.omega3_dri.day[i] <- b * o3
  reef_points$lost.zinc_dri.day[i] <- b * zn
  
}

# Capitalize countries
reef_points$country <- str_to_title(reef_points$country)

# Order by prevalence of overfishing
reef_points$country <- fct_relevel(reef_points$country, "Reunion", "Mozambique", "Kenya",
  "Tanzania", "Mauritius", "Madagascar", "Seychelles", "Mayotte", "Comoros", "Maldives")

# Function to make data for plotting simple donuts where eez is the
# name of the EEZ as it appears in reef_points. NB: requires reef_points
donuts <- function(eez){
  
  # Select for this country
  x <- filter(reef_points, country == eez)
  
  # Relabel underfished cells as sustainable
  x$status <- str_replace_all(x$status, pattern = "Underfished", replacement = "Sustainable")
  
  # New data frame of status
  y <- data.frame(table(x$status))
  
  # Rename columns
  colnames(y) <- c("status", "count")
  
  # Add country column
  y$country <- eez
  
  # Add additional columns for donut plotting
  y$fraction <- y$count / sum(y$count)
  y$ymax <- cumsum(y$fraction)
  y$ymin <- c(0, head(y$ymax, n=-1))
  
  return(y)
  
}

# Temp data for plotting donuts for each eez - changing underfished to sustainable
maldives <- donuts("Maldives")
comoros <- donuts("Comoros")
mozambique <- donuts("Mozambique")
mayotte <- donuts("Mayotte")
madagascar <- donuts("Madagascar")
mauritius <- donuts("Mauritius")
reunion <- donuts("Reunion")
kenya <- donuts("Kenya")
tanzania <- donuts("Tanzania")
seychelles <- donuts("Seychelles")

# Bind rows of donuts
temp.donuts <- bind_rows(maldives, comoros, mozambique, mayotte, madagascar, mauritius, reunion, kenya, tanzania, seychelles)

# Get order of countries' prevalence of overfishing
temp.donuts2 <- filter(temp.donuts, status == "Overfished")

# Order countries as factor
temp.donuts$country <- fct_relevel(temp.donuts$country, "Reunion", "Mozambique", "Kenya",
  "Tanzania", "Mauritius", "Madagascar", "Seychelles", "Mayotte", "Comoros", "Maldives")

# Plot donuts
ggplot(temp.donuts, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill = status)) +
  geom_rect() +
  facet_wrap("country", ncol = 3) +
  coord_polar(theta="y") + 
  xlim(c(2.3, 4)) +
  scale_fill_discrete(type = c("red", "blue")) +
  guides(fill = guide_legend(title = "", ncol = 1, label.position = "bottom",
    keywidth = unit(1, "mm"), keyheight = unit(4, "mm"))) +
  theme_void() +
  #labs(caption = "Prevalence of Overfished\nReef Sites") +
  theme(text = element_text(size = 6.5),
    plot.caption = element_text(hjust = 0.5),
    strip.text = element_text(size = 6.5),
    legend.position = c(0.85, 0.125),
    plot.margin = unit(c(5, 1, 1, 1), units = "mm")
    )

# Save image
ggsave("figures/05_PrevalenceOverfishing.pdf",
  width = 85, height = 100, units = "mm")


## Plot lost nutrient yields

# Subset data for overfished points only
temp <- filter(reef_points, status == "Overfished")

# Order countries as factor
temp$country <- fct_reorder(temp$country, temp$lost.yield_kg.day)

# Plot lost selenium yields
fig4a <- ggplot(temp, aes(x = lost.selenium_dri.day, y = country)) +
  geom_violin(fill = "gray", alpha = 0.5, color = "#0073C2FF",
    draw_quantiles = 0.5, position = "identity") +
  labs(title = "Selenium", color = "", y = "",
    x = "") +
  coord_cartesian(xlim = c(0, 100)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(legend.position = "none",
    text = element_text(size = 6.5),
    title = element_text(face = "bold"),
    plot.margin = unit(c(1, 3, 0, 0), units = "mm")
    )

# Plot lost omega-3 yields
fig4b <- ggplot(temp, aes(x = lost.omega3_dri.day, y = country)) +
  geom_violin(fill = "gray", alpha = 0.5, color = "#7AA6DCFF",
    draw_quantiles = 0.5, position = "identity") +
  labs(title = "Omega-3", color = "", y = "",
    x = "") +
  coord_cartesian(xlim = c(0, 100)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(legend.position = "none",
    text = element_text(size = 6.5),
    title = element_text(face = "bold"),
    plot.margin = unit(c(0, 3, 0, 0), units = "mm")
    )

# Plot lost zinc yields
fig4c <- ggplot(temp, aes(x = lost.zinc_dri.day, y = country)) +
  geom_violin(fill = "gray", alpha = 0.5, color = "#CD534CFF",
    draw_quantiles = 0.5, position = "identity") +
  labs(title = "Zinc", color = "", y = "",
    x = expression(paste("Lost Yield at Overfished Sites ", bgroup("(", frac('DRI', 'km'^2~day), ")"), sep = ""))) +
  coord_cartesian(xlim = c(0, 100)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(legend.position = "none",
    text = element_text(size = 6.5),
    title = element_text(face = "bold"),
    plot.margin = unit(c(0, 3, 0, 0), units = "mm")
    )

# Combine plots
ggarrange(fig4a, fig4b, fig4c,
  font.label = list(size = 6.5),
  nrow = 3, ncol = 1)

# Save figure
ggsave("figures/05_NutrientYieldsLost.pdf",
  width = 85, #height = 70,
  units = "mm",
  dpi = 300)



