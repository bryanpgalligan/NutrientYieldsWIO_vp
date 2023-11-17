## This script cleans and explores the detailed catch monitoring data and generates:
##    (1) national-level landings by species, averaged across sites (one file per country)
##    (2) WCS catch monitoring data aggregated to trip level to compare nutrition
##        outcomes to management classification




##### Load Data #####

# Species key downloaded from FishBase
species_data <- read.csv("data/clean-data/02_SpeciesData.csv")

# Dietary reference intakes
dri <- read.csv("data/clean-data/02_DietaryReferenceIntakes.csv")




##### Functions #####

## A function that pulls the WoRMS results you actually want: only the best matches
## for a list of taxa, returned as a dataframe with categories including Aphia ID,
## accepted name, supplied name, and full taxonomic information. This function can be
## used for lists of up to 350 species. NB: This is a wrapper for the
## wm_records_taxamatch() function in the worrms package. The wrapper will also load
## the necessary packages worrms and dplyr.

## Parameter: a list of species' scientific names

## Output: a dataframe of WoRMS results

worms_taxa2 <- function(species){
  
  # Load packages
  library(worrms)
  library(dplyr)
  
  ## Retrieve WoRMS results based on length of species list
  
  if (length(species) < 51){
    
    temp <- wm_records_taxamatch(species)
    df <- bind_rows(temp)
    
  }
  
  if (length(species) > 50 && length(species) < 101){
    
    temp <- wm_records_taxamatch(species[1:50])
    temp2 <- wm_records_taxamatch(species[51:length(species)])
    df <- bind_rows(temp, temp2)
    
  }
  
  if (length(species) > 100 && length(species) < 151){
    
    temp <- wm_records_taxamatch(species[1:50])
    temp2 <- wm_records_taxamatch(species[51:100])
    temp3 <- wm_records_taxamatch(species[101:length(species)])
    df <- bind_rows(temp, temp2, temp3)
    
  }
  
  if (length(species) > 150 && length(species) < 201){
    
    temp <- wm_records_taxamatch(species[1:50])
    temp2 <- wm_records_taxamatch(species[51:100])
    temp3 <- wm_records_taxamatch(species[101:150])
    temp4 <- wm_records_taxamatch(species[151:length(species)])
    df <- bind_rows(temp, temp2, temp3, temp4)
    
  }
  
  if (length(species) > 200 && length(species) < 251){
    
    temp <- wm_records_taxamatch(species[1:50])
    temp2 <- wm_records_taxamatch(species[51:100])
    temp3 <- wm_records_taxamatch(species[101:150])
    temp4 <- wm_records_taxamatch(species[151:200])
    temp5 <- wm_records_taxamatch(species[201:length(species)])
    df <- bind_rows(temp, temp2, temp3, temp4, temp5)
    
  }
  
  if (length(species) > 250 && length(species) < 301){
    
    temp <- wm_records_taxamatch(species[1:50])
    temp2 <- wm_records_taxamatch(species[51:100])
    temp3 <- wm_records_taxamatch(species[101:150])
    temp4 <- wm_records_taxamatch(species[151:200])
    temp5 <- wm_records_taxamatch(species[201:250])
    temp6 <- wm_records_taxamatch(species[251:length(species)])
    df <- bind_rows(temp, temp2, temp3, temp4, temp5, temp6)
    
  }
  
  if (length(species) > 300 && length(species) < 351){
    temp <- wm_records_taxamatch(species[1:50])
    temp2 <- wm_records_taxamatch(species[51:100])
    temp3 <- wm_records_taxamatch(species[101:150])
    temp4 <- wm_records_taxamatch(species[151:200])
    temp5 <- wm_records_taxamatch(species[201:250])
    temp6 <- wm_records_taxamatch(species[251:300])
    temp7 <- wm_records_taxamatch(species[301:length(species)])
    df <- bind_rows(temp, temp2, temp3, temp4, temp5, temp6, temp7)
  }
  
  ## Remove duplicates from df
  if (nrow(df) > 1){
    
    # Empty vector to hold list of row numbers to be deleted
    temp <- c()
    
    # Retrieve row numbers to be deleted
    for(i in 2:nrow(df)){
        
      # test if this provided name matches the one in the row above it
      if (df$scientificname[i] == df$scientificname[i-1]){
        
        # if it is a duplicate, add to temp vector
        temp <- append(temp, i, after = length(temp))
      }
        
    }
    
    # Ensure there were actually duplicates before deleting
    if(length(temp) > 0){

      # Delete duplicate rows from WoRMS results
      df <- df[-temp, ]
    
    }
    
  }
  
  
  # Select relevant columns from WoRMS results
  df <- select(df, valid_AphiaID, scientificname, valid_name, rank,
    kingdom, phylum, class, order, family, genus)
    
  # Rename columns as desired
  df <- rename(df, c(AphiaID = valid_AphiaID,
    supplied_name = scientificname))
    
  return(df)

}


## A function that generates random alphanumeric strings
randos <- function(n = 100000) {
  b <- do.call(paste0, replicate(5, sample(LETTERS, n, TRUE), FALSE))
  paste0(b, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
}

# Master list of trip id's
trips <- randos()


## A function that calculates fish weights based on measured length

## Parameters:
##    x:      a data frame of catch data that at least has the following columns:
##              species (can actually be a species, genus, or family)
##              length_cm (total length in cm)
##              weight_kg (fish weight in kg, can be empty)
##    key:    a data frame created by the estimate() function in rfishbase that at
##            least has columns:
##              species
##              a
##              b
##    worms:  results of the worms_taxa2 function

## Output: a data frame with an additional column for fish weight

## NB: if there are existing fish weights in the weight_kg column, they will be
## replaced if they can be estimated from LWR. If they can't be estimated, they will
## not be replaced. If you wish to preserve measured weights, you can save them in a
## new column under a different name. If you don't want to mix estimated and measured
## weights, you can delete all measured weights in the weight_kg column.

fish_weight <- function(x, key, worms){
    
  for (i in 1:nrow(x)){
    
    # Extract L
    L <- x$length_cm[i]
    
    # Find row number for this species in worms
    j <- match(x$species[i], worms$valid_name)
    
    # If ID is to species level
    if (worms$rank[j] == "Species"){
      
      # Find row number for this species in key
      k <- which(key$species %in% x$species[i])
      
      # And if this species has LWR estimates
      if (is.na(key$a[k]) == FALSE){
        
        # Extract a and b
        a <- key$a[k]
        b <- key$b[k]
        
        # Calculate weight based on the formula W = a * L^b and convert g to kg
        x$weight_kg[i] <- (a * L^b) / 1000
        
      } else {
        
        # If species does not have estimates for a and b
        
        # Subset the key for all species in this genus
        genus <- key[str_which(key$species, pattern = str_split_1(worms$valid_name[j], pattern = " ")[1]), ]
        
        # Extract mean values for a and b
        a <- mean(genus$a, na.rm = TRUE)
        b <- mean(genus$b, na.rm = TRUE)
        
        # If the genus has LWR estimates
        if (is.nan(a) == FALSE){
          
          # Calculate weight based on the formula W = a * L^b and convert g to kg
          x$weight_kg[i] <- (a * L^b) / 1000
          
        } else {
          
          # If the genus does not have estimates for a and b
          
          # Subset the key for all species in this family
          family <- filter(key, key$family == worms$family[j])
          
          # Extract mean values for a and b
          a <- mean(family$a, na.rm = TRUE)
          b <- mean(family$b, na.rm = TRUE)
          
          # Calculate weight based on the formula W = a * L^b and convert g to kg
          x$weight_kg[i] <- (a * L^b) / 1000
          
        }
        
      }
      
    }
    
    # If ID is to genus level
    if (worms$rank[j] == "Genus"){
        
      # Subset the key for all species in this genus
      genus <- key[str_which(key$species, pattern = str_split_1(worms$valid_name[j], pattern = " ")[1]), ]
        
      # Extract mean values for a and b
      a <- mean(genus$a, na.rm = TRUE)
      b <- mean(genus$b, na.rm = TRUE)
        
      # If the genus has LWR estimates
      if (is.nan(a) == FALSE){
        
        # Calculate weight based on the formula W = a * L^b and convert g to kg
        x$weight_kg[i] <- (a * L^b) / 1000
          
      } else {
          
        # If the genus does not have estimates for a and b
          
        # Subset the key for all species in this family
        family <- filter(key, key$family == worms$family[j])
          
        # Extract mean values for a and b
        a <- mean(family$a, na.rm = TRUE)
        b <- mean(family$b, na.rm = TRUE)
          
        # Calculate weight based on the formula W = a * L^b and convert g to kg
        x$weight_kg[i] <- (a * L^b) / 1000
      
      }
      
    }
    
    # If ID is to family level
    if(worms$rank[j] == "Family"){
      
        # Subset the key for all species in this family
        family <- filter(key, key$family == worms$family[j])
          
        # Extract mean values for a and b
        a <- mean(family$a, na.rm = TRUE)
        b <- mean(family$b, na.rm = TRUE)
          
        # Calculate weight based on the formula W = a * L^b and convert g to kg
        x$weight_kg[i] <- (a * L^b) / 1000
      
    }
    
  }
  
  return(x)
  
}




## A function that calculates nutrient weights based on fish weights
## and also adds K, trophic level, and habitat

## Parameters:
##    x:      a data frame that at least has the following columns:
##              species (can actually be a species, genus, or family)
##              weight_kg (fish weight in kg)
##              length_cm (fish length in cm)
##    key:    a data frame created by the table functions in rfishbase that at
##            least has columns for species name, nutrient concentrations,
##            trophic level, K (K can be sourced from FishBase or from
##            FishLife, as it is here), and habitat
##    worms:  results of the worms_taxa2 function

## Output: a data frame with additional columns for nutrient weights and L/Lopt

nutrient_weights <- function(x, key, worms){
  
  # Add new columns to x
  x$calcium_mg <- NA
  x$iron_mg <- NA
  x$omega.3_g <- NA
  x$selenium_ug <- NA
  x$vitamin.a_ug <- NA
  x$zinc_mg <- NA
  x$k <- NA
  x$troph <- NA
  x$habitat <- NA
  
  # Add new columns to worms
  worms$calcium_mg <- NA
  worms$iron_mg <- NA
  worms$omega.3_g <- NA
  worms$selenium_ug <- NA
  worms$vitamin.a_ug <- NA
  worms$zinc_mg <- NA
  worms$nut.rank <- worms$rank
  worms$k <- NA
  worms$troph <- NA
  worms$habitat <- NA

  # Populate nutrient columns of worms
  for (i in 1:nrow(worms)){
    
    # If taxon is a species
    if(worms$rank[i] == "Species"){
      
      # Does this species have nutrient estimates?
      if(is.na(key$calcium_mg.100g[key$species == worms$valid_name[i]]) == FALSE){
        
        # If yes...
        
        # Populate columns in worms
        worms$calcium_mg[i] <- key$calcium_mg.100g[key$species == worms$valid_name[i]]
        worms$iron_mg[i] <- key$iron_mg.100g[key$species == worms$valid_name[i]]
        worms$omega.3_g[i] <- key$omega3_g.100g[key$species == worms$valid_name[i]]
        worms$selenium_ug[i] <- key$selenium_ug.100g[key$species == worms$valid_name[i]]
        worms$vitamin.a_ug[i] <- key$vitamin.a_ug.100g[key$species == worms$valid_name[i]]
        worms$zinc_mg[i] <- key$zinc_mg.100g[key$species == worms$valid_name[i]]
        
      } else{ # If not...
        
        # Replace rank with Genus
        worms$nut.rank[i] <- "Genus"
        
      }
      
    }
    
    # If taxon requires a genus-level estimate
    if(worms$nut.rank[i] == "Genus"){
      
      # Subset key for this genus
      key2 <- key[str_which(key$species, pattern = str_split_1(worms$valid_name[i], pattern = " ")[1]), ]
      
      # If the genus has estimates
      if(is.nan(mean(key2$calcium_mg.100g, na.rm = TRUE)) == FALSE){
        
        # Populate columns in worms
        worms$calcium_mg[i] <- mean(key2$calcium_mg.100g, na.rm = TRUE)
        worms$iron_mg[i] <- mean(key2$iron_mg.100g, na.rm = TRUE)
        worms$omega.3_g[i] <- mean(key2$omega3_g.100g, na.rm = TRUE)
        worms$selenium_ug[i] <- mean(key2$selenium_ug.100g, na.rm = TRUE)
        worms$vitamin.a_ug[i] <- mean(key2$vitamin.a_ug.100g, na.rm = TRUE)
        worms$zinc_mg[i] <- mean(key2$zinc_mg.100g, na.rm = TRUE)
        
      } else{
        
        # Require a family level estimate
        worms$nut.rank[i] <- "Family"
        
      }
      
    }
      
      # If a family estimate is required
      if(worms$nut.rank[i] == "Family"){
        
        # Subset key for all species in family
        key2 <- filter(key, key$family == worms$family[i])
        
        # Populate columns in worms
        worms$calcium_mg[i] <- mean(key2$calcium_mg.100g, na.rm = TRUE)
        worms$iron_mg[i] <- mean(key2$iron_mg.100g, na.rm = TRUE)
        worms$omega.3_g[i] <- mean(key2$omega3_g.100g, na.rm = TRUE)
        worms$selenium_ug[i] <- mean(key2$selenium_ug.100g, na.rm = TRUE)
        worms$vitamin.a_ug[i] <- mean(key2$vitamin.a_ug.100g, na.rm = TRUE)
        worms$zinc_mg[i] <- mean(key2$zinc_mg.100g, na.rm = TRUE)
        
      }
    
    
    # Populate nutrient columns in x
    x$calcium_mg[x$species == worms$valid_name[i]] <- (worms$calcium_mg[i] * 10) * x$weight_kg[x$species == worms$valid_name[i]]
    x$iron_mg[x$species == worms$valid_name[i]] <- (worms$iron_mg[i] * 10) * x$weight_kg[x$species == worms$valid_name[i]]
    x$omega.3_g[x$species == worms$valid_name[i]] <- (worms$omega.3_g[i] * 10) * x$weight_kg[x$species == worms$valid_name[i]]
    x$selenium_ug[x$species == worms$valid_name[i]] <- (worms$selenium_ug[i] * 10) * x$weight_kg[x$species == worms$valid_name[i]]
    x$vitamin.a_ug[x$species == worms$valid_name[i]] <- (worms$vitamin.a_ug[i] * 10) * x$weight_kg[x$species == worms$valid_name[i]]
    x$zinc_mg[x$species == worms$valid_name[i]] <- (worms$zinc_mg[i] * 10) * x$weight_kg[x$species == worms$valid_name[i]]
    
    
    # Populate trophic level, habitat, and k columns of worms (there are no NAs in species_data)
    
    # If taxon is a species
    if(worms$rank[i] == "Species"){
        
      # Populate columns
      worms$troph[i] <- key$troph[key$species == worms$valid_name[i]]
      worms$habitat[i] <- key$habitat[key$species == worms$valid_name[i]]
      worms$k[i] <- key$k[key$species == worms$valid_name[i]]
        
    }
      
    # If taxon is a genus
    if(worms$rank[i] == "Genus"){
        
      # Subset key for this genus
      key2 <- key[str_which(key$species, pattern = str_split_1(worms$valid_name[i], pattern = " ")[1]), ]
        
      # Populate columns
      worms$troph[i] <- mean(key2$troph, na.rm = TRUE)
      worms$habitat[i] <- names(which.max(table(key2$habitat)))
      worms$k[i] <- mean(key2$k, na.rm = TRUE)
        
    }
      
    # If taxon is a family
    if(worms$rank[i] == "Family"){
        
      # Subset key for all species in family
      key2 <- filter(key, key$family == worms$family[i])
        
      # Populate columns
      worms$troph[i] <- mean(key2$troph, na.rm = TRUE)
      worms$habitat[i] <- names(which.max(table(key2$habitat)))
      worms$k[i] <- mean(key2$k, na.rm = TRUE)
        
    }
    
    
    # Populate life history columns in x
    x$troph[x$species == worms$valid_name[i]] <- worms$troph[i]
    x$habitat[x$species == worms$valid_name[i]] <- worms$habitat[i]
    x$k[x$species == worms$valid_name[i]] <- worms$k[i]
    
  }
    
  return(x)
  
}




## A function that adds a column for the ratio of length to optimum length (L/Lopt)
## to detailed catch monitoring data.

## Parameters:
##      x:      a table that includes the following columns:
##                  species:    taxonomic identity of individual fish (either species,
##                              genus, or family) - scientific name
##                  length_cm:  length of fish in cm
##      worms:  results of the worms_taxa2 function for all taxa listed in x
##      key:    a table that includes the following columns:
##                  species:    scientific names of all species included in x as well
##                              as the constituent species of genuses and families in x
##                  Lopt_cm:    an estimate of optimum length (Lopt) in cm for each
##                              species. This is best derived from Jim Thorsen's
##                              FishLife package.

lopt_indicator <- function(x, key, worms){
  
  # Add additional columns
  x$LLopt <- NA
  x$LLmat <- NA
  worms$Lopt <- NA
  worms$Lmat <- NA
  worms$loptrank <- worms$rank
  
  for (i in 1:nrow(worms)){
    
    # If taxon is a species
    if(worms$rank[i] == "Species"){
      
      # Does this species have an Lopt estimate?
      if(is.na(key$Lopt_cm[key$species == worms$valid_name[i]]) == FALSE){
        
        # If so...
        
        # Populate worms
        worms$Lopt[i] <- key$Lopt_cm[key$species == worms$valid_name[i]]
        worms$Lmat[i] <- key$Lmat_cm[key$species == worms$valid_name[i]]
        
      } else{ # if not...
        
        # Change estimation rank to Genus
        worms$loptrank[i] <- "Genus"
        
      }
      
    }
    
    # If a genus-level estimate is needed
    if(worms$loptrank[i] == "Genus"){
      
      # Subset key for this genus
      key2 <- key[str_which(key$species, pattern = str_split_1(worms$valid_name[i], pattern = " ")[1]), ]
      
      # Populate worms
      worms$Lopt[i] <- mean(key2$Lopt_cm, na.rm = TRUE)
      worms$Lmat[i] <- mean(key2$Lmat_cm, na.rm = TRUE)
      
    }
    
    # If a family-level estimate is needed
    if(worms$loptrank[i] == "Family"){
      
      # Subset key for this family
      key2 <- filter(key, key$family == worms$family[i])
      
      # Populate worms
      worms$Lopt[i] <- mean(key2$Lopt_cm, na.rm = TRUE)
      worms$Lmat[i] <- mean(key2$Lmat_cm, na.rm = TRUE)
      
    }
    
    # Calculate L/Lopt
    x$LLopt[x$species == worms$valid_name[i]] <-
      x$length_cm[x$species == worms$valid_name[i]] / worms$Lopt[i]
    
    # Calculate L/Lmat
    x$LLmat[x$species == worms$valid_name[i]] <-
      x$length_cm[x$species == worms$valid_name[i]] / worms$Lmat[i]
    
  }
  
  return(x)
  
}




## A function that transforms detailed catch monitoring data (in which each row
## corresponds to one individual fish) into trip-level data.

## Parameters:
##      x: a table that includes the following columns:
##              trip:         unique identifier for each trip
##              source:       source for these data
##              date:         date of trip
##              country:      location of landing site
##              site:         name of landing site
##              effort:       number of crew, number of gears, or other measure of effort
##              gear:         fishing gear used
##              management:   management regime at this landing site
##              weight_kg:    weight of individual fish in kg
##              calcium_mg:   calcium in mg
##              iron_mg:      iron in mg
##              omega.3_g:    omega-3 in g
##              selenium_ug:  selenium in ug
##              vitamin.a_ug: vitamin a in ug
##              zinc_mg:      zinc in mg
##              LLopt:        ratio of L to Lopt
##      dri: daily recommended intakes for micronutrients with two columns
##              nutrients:  list of nutrients, including those listed in x
##              dri:        daily recommended intakes for each nutrient

## Output: a consolidated table with trip-level data

fish2trip <- function(x, dri){
  
  # Sum landings and nutrient weights by trip
  temp <- x %>%
    group_by(trip) %>%
    summarize(across(c(weight_kg, calcium_mg, iron_mg, omega.3_g, selenium_ug,
      vitamin.a_ug, zinc_mg), sum))
  
  # Rename weight to landings
  temp <- rename(temp, landings_kg = weight_kg)
  
  # Add additional columns
  temp$source <- NA
  temp$date <- NA
  temp$country <- NA
  temp$site <- NA
  temp$gear <- NA
  temp$management <- NA
  temp$effort <- NA
  temp$cpue_kg.effort.day <- NA
  temp$calcium_mg.pue <- NA
  temp$calcium_pdv.pue <- NA
  temp$iron_mg.pue <- NA
  temp$iron_pdv.pue <- NA
  temp$omega.3_g.pue <- NA
  temp$omega.3_pdv.pue <- NA
  temp$selenium_ug.pue <- NA
  temp$selenium_pdv.pue <- NA
  temp$vitamin.a_ug.pue <- NA
  temp$zinc_mg.pue <- NA
  temp$zinc_pdv.pue <- NA
  temp$species_richness <- NA
  temp$nutrient_evenness <- NA
  temp$LLmat <- NA
  temp$LLopt <- NA
  temp$length_cm <- NA
  temp$k <- NA
  temp$troph <- NA
  temp$habitat <- NA
  temp$nutrients_pdv.pue <- NA
  
  # Populate additional columns
  for (i in 1:nrow(temp)){
    
    # Extract index of the first row of this trip in x
    a <- match(temp$trip[i], x$trip)
    
    # Date
    temp$date[i] <- x$date[a]
    
    # Site
    temp$site[i] <- x$site[a]
    
    # Effort
    temp$effort[i] <- x$effort[a]
    
    # Gear
    temp$gear[i] <- x$gear[a]
    
    # Management
    temp$management[i] <- x$management[a]
    
    # Country
    temp$country[i] <- x$country[a]
    
    # Source
    temp$source[i] <- x$source[a]
    
  }
  
  # Calculate species richness, LLopt, LLmat, and length
  for (i in 1:nrow(temp)){
    
    # Subset of x for this trip
    y <- filter(x, trip == temp$trip[i])
    
    # Count unique species and divide by biomass
    z <- length(unique(y$species))
    z <- z / sum(y$weight_kg)
    
    # Save to temp
    temp$species_richness[i] <- z
    
    # Find mean LLopt
    z <- mean(y$LLopt, na.rm = TRUE)
    
    # Save to temp
    temp$LLopt[i] <- z
    
    # Find mean LLmat
    z <- mean(y$LLmat, na.rm = TRUE)
    
    # Save to temp
    temp$LLmat[i] <- z
    
    # Find mean length
    z <- mean(y$length_cm, na.rm = TRUE)
    
    # Save to temp
    temp$length_cm[i] <- z
    
    # Find mean k
    z <- mean(y$k, na.rm = TRUE)
    
    # Save to temp
    temp$k[i] <- z
    
    # Find mean troph
    z <- mean(y$troph, na.rm = TRUE)
    
    # Save to temp
    temp$troph[i] <- z
    
    # Most frequent habitat
    z <- names(which.max(table(y$habitat)))
    
    # Save to temp if there is anything to save
    if (length(z > 0)){
      temp$habitat[i] <- z
    }
    
  }
  
  # Calculate catch and nutrients per unit effort
  temp$cpue_kg.effort.day <- temp$landings_kg / temp$effort
  temp$calcium_mg.pue <- temp$calcium_mg / temp$effort
  temp$iron_mg.pue <- temp$iron_mg / temp$effort
  temp$omega.3_g.pue <- temp$omega.3_g / temp$effort
  temp$selenium_ug.pue <- temp$selenium_ug / temp$effort
  temp$vitamin.a_ug.pue <- temp$vitamin.a_ug / temp$effort
  temp$zinc_mg.pue <- temp$zinc_mg / temp$effort
  
  # Calculate nutrients percent daily value per unit effort
  temp$calcium_pdv.pue <- temp$calcium_mg.pue / dri[dri$nutrients == "calcium_mg", "dri"]
  temp$iron_pdv.pue <- temp$iron_mg.pue / dri[dri$nutrients == "iron_mg", "dri"]
  temp$omega.3_pdv.pue <- temp$omega.3_g.pue / dri[dri$nutrients == "omega3_g", "dri"]
  temp$selenium_pdv.pue <- temp$selenium_ug.pue / dri[dri$nutrients == "selenium_ug", "dri"]
  temp$vitamin.a_pdv.pue <- temp$vitamin.a_ug.pue / dri[dri$nutrients == "vitamina_ug", "dri"]
  temp$zinc_pdv.pue <- temp$zinc_mg.pue / dri[dri$nutrients == "zinc_mg", "dri"]
  
  
  ## Add columns for percent daily values and total nutrient density (out of 600%)

  # Calcium
  temp$calcium_mg.per.100g <- (temp$calcium_mg / temp$landings_kg) / 10
  temp$calcium_pdv <- temp$calcium_mg.per.100g / dri[dri$nutrients == "calcium_mg", "dri"]
  
  # Iron
  temp$iron_mg.per.100g <- (temp$iron_mg / temp$landings_kg) / 10
  temp$iron_pdv <- temp$iron_mg.per.100g / dri[dri$nutrients == "iron_mg", "dri"]
  
  # Omega 3
  temp$omega.3_g.per.100g <- (temp$omega.3_g / temp$landings_kg) / 10
  temp$omega.3_pdv <- temp$omega.3_g.per.100g / dri[dri$nutrients == "omega3_g", "dri"]
  
  # Selenium
  temp$selenium_ug.per.100g <- (temp$selenium_ug / temp$landings_kg) / 10
  temp$selenium_pdv <- temp$selenium_ug.per.100g / dri[dri$nutrients == "selenium_ug", "dri"]
  
  # Vitamin A
  temp$vitamin.a_ug.per.100g <- (temp$vitamin.a_ug / temp$landings_kg) / 10
  temp$vitamin.a_pdv <- temp$vitamin.a_ug.per.100g / dri[dri$nutrients == "vitamina_ug", "dri"]
  
  # Zinc
  temp$zinc_mg.per.100g <- (temp$zinc_mg / temp$landings_kg) / 10
  temp$zinc_pdv <- temp$zinc_mg.per.100g / dri[dri$nutrients == "zinc_mg", "dri"]
  
  # Cap percent daily values at 1
  for (i in 1:nrow(temp)){
    
    # If there are nutrients for this trip
    if(is.na(temp$calcium_pdv[i]) == FALSE){
    
      if (temp$calcium_pdv[i] > 1){
        temp$calcium_pdv[i] <- 1
      }
      
      if (temp$iron_pdv[i] > 1){
        temp$iron_pdv[i] <- 1
      }
      
      if (temp$omega.3_pdv[i] > 1){
        temp$omega.3_pdv[i] <- 1
      }
      
      if (temp$selenium_pdv[i] > 1){
        temp$selenium_pdv[i] <- 1
      }
      
      if (temp$vitamin.a_pdv[i] > 1){
        temp$vitamin.a_pdv[i] <- 1
      }
      
      if (temp$zinc_pdv[i] > 1){
        temp$zinc_pdv[i] <- 1
      }
      
    }
    
  }
  
  # Total nutrient composition out of 6 (600%)
  temp$nutrients_pdv <- temp$calcium_pdv + temp$iron_pdv +
    temp$omega.3_pdv + temp$selenium_pdv + temp$vitamin.a_pdv +
    temp$zinc_pdv
  
  # Total nutrient yields (number of people's nutritional requirements met per fisher per day)
  temp$nutrients_pdv.pue <- ((temp$nutrients_pdv * 10) * temp$cpue_kg.effort.day) / 6
  
  
  # Calculate nutrient evenness
  for (i in 1:nrow(temp)){
    
    # Subset of x for this trip
    y <- filter(x, trip == temp$trip[i])
    
    # If there are nutrient values, continue
    if(is.nan(mean(y$calcium_mg, na.rm = TRUE)) == FALSE){
    
      # Count unique species
      z <- length(unique(y$species))
      
      ## Add columns for percent daily values and total nutrients relative to dri
    
      # Calcium
      y$calcium_pdv <- y$calcium_mg / dri[dri$nutrients == "calcium_mg", "dri"]
      
      # Iron
      y$iron_pdv <- y$iron_mg / dri[dri$nutrients == "iron_mg", "dri"]
      
      # Omega 3
      y$omega.3_pdv <- y$omega.3_g / dri[dri$nutrients == "omega3_g", "dri"]
      
      # Selenium
      y$selenium_pdv <- y$selenium_ug / dri[dri$nutrients == "selenium_ug", "dri"]
      
      # Vitamin A
      y$vitamin.a_pdv <- y$vitamin.a_ug / dri[dri$nutrients == "vitamina_ug", "dri"]
      
      # Zinc
      y$zinc_pdv <- y$zinc_mg / dri[dri$nutrients == "zinc_mg", "dri"]
      
      # Total nutrient composition out of 6 (600%)
      y$nutrients_pdv <- y$calcium_pdv + y$iron_pdv +
        y$omega.3_pdv + y$selenium_pdv + y$vitamin.a_pdv +
        y$zinc_pdv
      
      
      ## Calculate nutrient evenness
      
      # Aggregate y by nutrients
      y <- aggregate(nutrients_pdv ~ species, y, FUN = sum)
      
      # Proportion (p)
      y$p <- y$nutrients_pdv / sum(y$nutrients_pdv)
      
      # Natural log of proportion
      y$lnp <- log(y$p)
      
      # p * ln(p)
      y$plnp <- y$p * y$lnp
      
      # Shannon diversity index (h)
      h <- -1 * sum(y$plnp)
      
      # Calculate evenness
      e <- h / log(z)
      
      # Save to temp
      temp$nutrient_evenness[i] <- e
   
    }
       
  }
  
  return(temp)
  
}




##### Dzoga 2020 #####

## This section cleans data originally used for:
## Dzoga, M., Simatele, D., & Munga, C. (2020). Characterisation of artisanal catches in
##    selected fishing areas of the Lower Tana Delta and Malindi-Ungwana Bay, Kenya. 
##    Western Indian Ocean Journal of Marine Science, 19(1), Article 1. 
##    https://doi.org/10.4314/wiojms.v19i1.4

# Import catch data
catch_kenya_dzoga <- read_excel("data/raw-data/catch-monitoring/kenya-catch-monitoring-dzoga2020.xlsx", 
    col_types = c("text", "text", "text", 
        "text", "text", "text", "text", "text", 
        "numeric", "text", "text", "numeric", 
        "text", "text", "numeric", "numeric", 
        "numeric", "skip"))


# Select desired columns
catch_kenya_dzoga <- select(catch_kenya_dzoga, "Date", "Crew no.", "Station",
  "Gear used", "Species name", "Sample wgt. In kgs.", "Lengths in cm.")

# Correct column names
catch_kenya_dzoga <- rename(catch_kenya_dzoga,
  c(date = "Date",
    effort = "Crew no.",
    site = "Station",
    gear = "Gear used",
    species = "Species name",
    weight_kg = "Sample wgt. In kgs.",
    length_cm = "Lengths in cm."))


## Site column

# Make lower case
catch_kenya_dzoga$site <- str_to_lower(catch_kenya_dzoga$site)

# Delete freshwater site
catch_kenya_dzoga <- filter(catch_kenya_dzoga, site != "ozi")


## Correct species names

# Change local names to scientific names
catch_kenya_dzoga$species <- str_replace_all(catch_kenya_dzoga$species,
  "Kidara", "Carangoides equula")
catch_kenya_dzoga$species <- str_replace_all(catch_kenya_dzoga$species,
  "Songoro", "Rachycentron canadum")
catch_kenya_dzoga$species <- str_replace_all(catch_kenya_dzoga$species,
  "Borode", "Chanos chanos")

# Delete freshwater species
catch_kenya_dzoga <- filter(catch_kenya_dzoga, species != "Clarias gariepinus")

# Make all scientific names sentence case
catch_kenya_dzoga$species <- str_to_sentence(catch_kenya_dzoga$species)

# Remove white space
catch_kenya_dzoga$species <- str_squish(catch_kenya_dzoga$species)

# Find missepllings with taxize package
species_dzoga <- unique(catch_kenya_dzoga$species)
temp <- gnr_resolve(species_dzoga, best_match_only = TRUE, canonical = TRUE)
temp <- filter(temp, temp$user_supplied_name != temp$matched_name2)

# Replace misspelled names
catch_kenya_dzoga$species <- str_replace_all(catch_kenya_dzoga$species,
  "Uroteuthis duvaucellii", "Uroteuthis duvaucelii")
catch_kenya_dzoga$species <- str_replace_all(catch_kenya_dzoga$species,
  "Carangoides fedau", "Carangoides ferdau")
catch_kenya_dzoga$species <- str_replace_all(catch_kenya_dzoga$species,
  "Cynoglosus lida", "Cynoglossus lida")
catch_kenya_dzoga$species <- str_replace_all(catch_kenya_dzoga$species,
  "Scylla cerrata", "Scylla serrata")
catch_kenya_dzoga$species <- str_replace_all(catch_kenya_dzoga$species,
  "Terapon jabuar", "Terapon jarbua")
catch_kenya_dzoga$species <- str_replace_all(catch_kenya_dzoga$species,
  "Liza vaigeinsis", "Liza vaigensis")
catch_kenya_dzoga$species <- str_replace_all(catch_kenya_dzoga$species,
  "Pelona ditchela", "Pellona ditchela")
catch_kenya_dzoga$species <- str_replace_all(catch_kenya_dzoga$species,
  "Trachinotus tol", "Trachinotus")
catch_kenya_dzoga$species <- str_replace_all(catch_kenya_dzoga$species,
  "Synoglosus lida", "Cynoglossus lida")
catch_kenya_dzoga$species <- str_replace_all(catch_kenya_dzoga$species,
  "Plectorhincus", "Plectorhinchus")
catch_kenya_dzoga$species <- str_replace_all(catch_kenya_dzoga$species,
  "Lethrinus rivulatus", "Lutjanus rivulatus")
catch_kenya_dzoga$species <- str_replace_all(catch_kenya_dzoga$species,
  "Lemipterous randali", "Nemipterus randalli")
catch_kenya_dzoga$species <- str_replace_all(catch_kenya_dzoga$species,
  "Trichiurus lepterus", "Trichiurus lepturus")

# Delete unresolvable names
catch_kenya_dzoga <- filter(catch_kenya_dzoga, species != "Nyesa")
catch_kenya_dzoga <- filter(catch_kenya_dzoga, species != "Sharui")
catch_kenya_dzoga <- filter(catch_kenya_dzoga, species != "Mbekwele")

# Get WoRMS matches
species_dzoga <- unique(catch_kenya_dzoga$species)
worms_dzoga <- worms_taxa2(species_dzoga)

# Correct species names in catch data
for (i in 1:nrow(catch_kenya_dzoga)){
  
  # Get index of row in WoRMS results
  x <- which(worms_dzoga$supplied_name == catch_kenya_dzoga$species[i])
  
  # Replace supplied name with new name
  catch_kenya_dzoga$species[i] <- worms_dzoga$valid_name[x]
  
}


# Remove non-chordata
temp <- subset(worms_dzoga$valid_name, worms_dzoga$phylum != "Chordata")
catch_kenya_dzoga <- filter(catch_kenya_dzoga, !(catch_kenya_dzoga$species %in% temp))
worms_dzoga <- filter(worms_dzoga, phylum == "Chordata")


## Add unique trip id's to each trip

# Empty column for trip ID
catch_kenya_dzoga$trip <- NA

# Fill in trip column
for (i in 1:nrow(catch_kenya_dzoga)){

  # If there is a crew size listed in this row
  if (is.na(catch_kenya_dzoga$effort[i]) == FALSE){

    # Add a trip ID
    catch_kenya_dzoga$trip[i] <- trips[i]

    # Remove trip ID from master list
    trips <- trips[-i]

  }

  # If there is no crew size listed in this row
  if (is.na(catch_kenya_dzoga$effort[i]) == TRUE){

    # Paste the previous trip ID
    catch_kenya_dzoga$trip[i] <- catch_kenya_dzoga$trip[i-1]

  }

}


## Date column

# Replace . with /
catch_kenya_dzoga$date <- str_replace_all(catch_kenya_dzoga$date, "\\.", "/")

# Temporary vector to hold new dates
temp <- c()

# Reformat dates
for (i in 1:nrow(catch_kenya_dzoga)){

  # If the date is in excel format
  if (str_length(catch_kenya_dzoga$date[i]) == 5){

    # Convert to numeric
    x <- as.numeric(catch_kenya_dzoga$date[i])

    # Convert to date
    x <- as.Date(x, origin = "1899-12-30")

    # Save to the temporary vector
    temp <- append(temp, x, after = length(temp))

  } else{

    # If the date is in some dd/mm/yyyy format
    x <- as.Date(catch_kenya_dzoga$date[i], format = "%d/%m/%Y")

    # Save to the temporary vector
    temp <- append(temp, x, after = length(temp))

  }

}

# Replace old date column with new one
catch_kenya_dzoga$date <- temp

# Delete all dates later than 2017 (these are 155 typos)
catch_kenya_dzoga <- filter(catch_kenya_dzoga, date < "2018-01-01")


## Country column
catch_kenya_dzoga$country <- "kenya"


## Source column
catch_kenya_dzoga$source <- "dzoga2020"


## Management column
catch_kenya_dzoga$management <- NA


## Gear column

# Lower case
catch_kenya_dzoga$gear <- str_to_lower(catch_kenya_dzoga$gear)

# Correct entries
catch_kenya_dzoga$gear <- str_replace_all(catch_kenya_dzoga$gear,
  pattern = "basket trap",
  replacement = "trap")
catch_kenya_dzoga$gear <- str_replace_all(catch_kenya_dzoga$gear,
  pattern = "monofilament net",
  replacement = "monofilament")
catch_kenya_dzoga$gear <- str_replace_all(catch_kenya_dzoga$gear,
  pattern = "seine net",
  replacement = "beachseine")
catch_kenya_dzoga$gear <- str_replace_all(catch_kenya_dzoga$gear,
  pattern = "speargun",
  replacement = "spear")
catch_kenya_dzoga$gear <- str_replace_all(catch_kenya_dzoga$gear,
  pattern = "ring net",
  replacement = "ringnet")


## Weight column

# Unfortunately, it looks like many of Dzoga's weights are not possible values. With
# no way to distinguish what is real and what is fake, I will replace all weights
# using published length-weight relationships obtained from FishBase.
catch_kenya_dzoga <- fish_weight(catch_kenya_dzoga, key = species_data,
  worms = worms_dzoga)


## Length column

# Delete rows with impossibly long fish (after manual inspection of the data)
for (i in 1:nrow(catch_kenya_dzoga)){
  
  if (catch_kenya_dzoga$species[i] == "Argyrops spinifer"){
    
    if (catch_kenya_dzoga$length_cm[i] > 70){
      
      catch_kenya_dzoga$length_cm[i] <- NA
      
    }
    
  }
  
  if (catch_kenya_dzoga$species[i] == "Caesio lunaris"){
    
    if (catch_kenya_dzoga$length_cm[i] > 40){
      
      catch_kenya_dzoga$length_cm[i] <- NA
      
    }
    
  }
  
}

# Delete rows with 0 or NA for length (and thus also weight)
catch_kenya_dzoga <- filter(catch_kenya_dzoga, length_cm > 0)


## Add Dzoga 2020 to Kenya catch summary

# Aggregate by landing site and species
temp <- aggregate(weight_kg ~ species + site,
  data = catch_kenya_dzoga, FUN = sum)

# Make a summary of landings by species by site
temp <- temp %>%
  pivot_wider(names_from = site, values_from = weight_kg)

# Save catch_kenya
catch_kenya <- temp


## Add Dzoga 2020 to trip data

# Calculate nutrient weights
catch_kenya_dzoga <- nutrient_weights(catch_kenya_dzoga, species_data, worms_dzoga)

# Add L/Lopt indicator
catch_kenya_dzoga <- lopt_indicator(catch_kenya_dzoga, species_data, worms_dzoga)

# Calculate trip data
trip_dzoga <- fish2trip(catch_kenya_dzoga, dri = dri)

# Select and reorder columns for master trip data
trip <- trip_dzoga %>% select(trip, date, country, site, management, gear,
  landings_kg, cpue_kg.effort.day, LLmat, LLopt, length_cm,
  k, troph, habitat, species_richness,
  nutrients_pdv, nutrients_pdv.pue, nutrient_evenness,
  calcium_mg, calcium_mg.pue, calcium_mg.per.100g, calcium_pdv, calcium_pdv.pue,
  iron_mg, iron_mg.pue, iron_mg.per.100g, iron_pdv, iron_pdv.pue,
  omega.3_g, omega.3_g.pue, omega.3_g.per.100g, omega.3_pdv, omega.3_pdv.pue,
  selenium_ug, selenium_ug.pue, selenium_ug.per.100g, selenium_pdv, selenium_pdv.pue,
  vitamin.a_ug, vitamin.a_ug.pue, vitamin.a_ug.per.100g, vitamin.a_pdv, vitamin.a_pdv.pue,
  zinc_mg, zinc_mg.pue, zinc_mg.per.100g, zinc_pdv, zinc_pdv.pue,
  source)

# Make date column in date format
trip$date <- as.Date(trip$date, origin = "1970-01-01")


## Save to fish master
fish <- catch_kenya_dzoga %>% select(trip, date, country, site, gear, management,
  species, weight_kg, length_cm, LLopt, LLmat, k, troph, habitat, source)



##### Musembi 2019 #####

## This section cleans the data originally used for:
## Musembi, P., Fulanda, B., Kairo, J., & Githaiga, M. (2019). Species composition, 
##    abundance, and fishing methods of small-scale fisheries in the seagrass meadows
##    of Gazi Bay, Kenya. Journal of the Indian Ocean Region, 15(2), 139–156. 
##    https://doi.org/10.1080/19480881.2019.1603608

## Musembi et al. used number of gears as a measure of effort instead of number
##    of fishers, so we will only use this for catch composition data.

# Load catch monitoring data
catch_kenya_musembi <- read_excel("data/raw-data/catch-monitoring/kenya-catch-monitoring-musembi2019.xlsx")

# Select only desired columns
catch_kenya_musembi <- select(catch_kenya_musembi, "Species", "Lengths (cm)",
  "Individual weight (kg)", "Catch event", "Date", "Gear type")

# Rename columns
catch_kenya_musembi <- rename(catch_kenya_musembi,
  species = "Species",
  length_cm = "Lengths (cm)",
  weight_kg = "Individual weight (kg)",
  trip = "Catch event",
  date = "Date",
  gear = "Gear type")


## Site column
catch_kenya_musembi$site <- "gazi"


## Species column

# Match taxa with worms
species_musembi <- unique(catch_kenya_musembi$species)
worms_musembi <- worms_taxa2(species_musembi)

# Replace incorrect species names
for (i in 1:nrow(worms_musembi)){
  
  if (worms_musembi$supplied_name[i] != worms_musembi$valid_name[i]){
    
    catch_kenya_musembi$species <- gsub(pattern = worms_musembi$supplied_name[i],
      replacement = worms_musembi$valid_name[i],
      catch_kenya_musembi$species,
      fixed = TRUE)
    
  }
  
}

# NB it failed to insert "Uroteuthis (Photololigo) duvaucelii" because of the
# parentheses, but we're going to delete that species anyway.

# Temporary data frame of non-chordata
temp <- filter(worms_musembi, phylum != "Chordata")
temp2 <- c()

# List of row numbers of non-chordata
for (i in 1:nrow(catch_kenya_musembi)){
  
  if (any(temp == catch_kenya_musembi$species[i]) == TRUE){
    
    temp2 <- append(temp2, i, after = length(temp2))
    
  }
  
}

# Delete non-chordata
catch_kenya_musembi <- catch_kenya_musembi[-temp2,]
worms_musembi <- worms_musembi[worms_musembi$phylum == "Chordata",]

## NB I have no idea why, but it still hasn't deleted Loligo duvaucelli. Let's delete
## that species manually.
catch_kenya_musembi <- filter(catch_kenya_musembi, species != "Loligo duvaucelli")

# Manually correct misspelled species
catch_kenya_musembi$species <- str_replace_all(catch_kenya_musembi$species,
  pattern = "Scarus psittatus", replacement = "Scarus psittacus")


## Trip column

# List of unique trip id's in existing data
oldtrips <- unique(catch_kenya_musembi$trip)

# Replace old trip id's with new ones
for (i in 1:length(oldtrips)){
  catch_kenya_musembi$trip[catch_kenya_musembi$trip == oldtrips[i]] <- trips[i]
}

# Delete used trip id's from trips vector
trips <- trips[-(1:length(oldtrips))]


## Date column
catch_kenya_musembi$date <- as.Date(catch_kenya_musembi$date)


## Gear column

# Make all lower case
catch_kenya_musembi$gear <- str_to_lower(catch_kenya_musembi$gear)

# Traps
catch_kenya_musembi$gear[catch_kenya_musembi$gear == "basket traps"] <- "trap"

# Cast nets
catch_kenya_musembi$gear[catch_kenya_musembi$gear == "cast net"] <- "castnet"

# Hook and stick
catch_kenya_musembi$gear[catch_kenya_musembi$gear == "hook and stick"] <- "hookandstick"

# Reef seine
catch_kenya_musembi$gear[catch_kenya_musembi$gear == "reef seine"] <- "reefseine"

# Spear gun
catch_kenya_musembi$gear[catch_kenya_musembi$gear == "spear gun"] <- "spear"


## Add nutrient weights and LLopt
catch_kenya_musembi <- nutrient_weights(catch_kenya_musembi, species_data, worms_musembi)
catch_kenya_musembi <- lopt_indicator(catch_kenya_musembi, species_data, worms_musembi)


## Generate and save trip data

# Add additional columns
catch_kenya_musembi$effort <- NA
catch_kenya_musembi$management <- NA
catch_kenya_musembi$country <- "kenya"
catch_kenya_musembi$source <- "musembi2019"

# Generate trip data
trip_musembi <- fish2trip(catch_kenya_musembi, dri = dri)

# Select and reorder columns
trip_musembi <- trip_musembi %>% select(trip, date, country, site, management, gear,
  landings_kg, cpue_kg.effort.day,
  LLmat, LLopt, length_cm,
  k, troph, habitat, species_richness,
  nutrients_pdv, nutrients_pdv.pue, nutrient_evenness,
  calcium_mg, calcium_mg.pue, calcium_mg.per.100g, calcium_pdv, calcium_pdv.pue,
  iron_mg, iron_mg.pue, iron_mg.per.100g, iron_pdv, iron_pdv.pue,
  omega.3_g, omega.3_g.pue, omega.3_g.per.100g, omega.3_pdv, omega.3_pdv.pue,
  selenium_ug, selenium_ug.pue, selenium_ug.per.100g, selenium_pdv, selenium_pdv.pue,
  vitamin.a_ug, vitamin.a_ug.pue, vitamin.a_ug.per.100g, vitamin.a_pdv, vitamin.a_pdv.pue,
  zinc_mg, zinc_mg.pue, zinc_mg.per.100g, zinc_pdv, zinc_pdv.pue,
  source)

# Make date column in date format
trip_musembi$date <- as.Date(trip_musembi$date, origin = "1970-01-01")

# Bind to master trip data
trip <- bind_rows(trip, trip_musembi)


## Add Musembi 2019 to catch_kenya

# Aggregate by landing site and species
temp <- aggregate(weight_kg ~ species + site,
  data = catch_kenya_musembi, FUN = sum)

# Make a summary of landings by species by site
temp <- temp %>%
  pivot_wider(names_from = site, values_from = weight_kg)

# Bind to catch_kenya master
catch_kenya <- full_join(catch_kenya, temp)


## Save to fish master
fish <- full_join(fish, catch_kenya_musembi)
fish <- select(fish, -(calcium_mg:effort))




##### WCS Kenya #####

## This section cleans the WCS catch monitoring data

# Load data
catch_kenya_wcs <- read_excel("data/raw-data/catch-monitoring/kenya-catch-monitoring-wcs.xlsx", 
    col_types = c("text", "numeric", "text", 
        "text", "numeric", "text", "text", 
        "text", "text", "text", "text", "text", 
        "text", "text", "text", "text", 
        "text", "text", "text", 
        "text", "text", "text", 
        "text", "text", "text", 
        "text", "text", "text", 
        "text", "text", "text", 
        "text", "text", "text", 
        "text", "text", "text", 
        "text", "text", "text", 
        "text", "text", "text", 
        "text", "text"))

# Rename columns
catch_kenya_wcs <- rename(catch_kenya_wcs,
  c(date = Date,
    year = Year,
    site = Landing,
    management = "Management type",
    crew_id = Group,
    crew_size = "#Fishers",
    percent_sampled = "% of catch sampled",
    gear = Fishgear,
    species = Name))


## Date column

# Temporary vector to hold new dates
temp <- c()

# Reformat dates
for (i in 1:nrow(catch_kenya_wcs)){
  
  # If the date is in excel format
  if (str_length(catch_kenya_wcs$date[i]) == 5){
    
    # Convert to numeric
    x <- as.numeric(catch_kenya_wcs$date[i])
    
    # Convert to date
    x <- as.Date(x, origin = "1899-12-30")
    
    # Save to the temporary vector
    temp <- append(temp, x, after = length(temp))
    
  } else{
    
    # If the date is in some dd/mm/yyyy format
    x <- as.Date(catch_kenya_wcs$date[i], format = "%d/%m/%Y")
    
    # Save to the temporary vector
    temp <- append(temp, x, after = length(temp))
    
  }
  
}

# Replace old date column with new one
catch_kenya_wcs$date <- temp


## Site column
catch_kenya_wcs$site <- str_to_lower(catch_kenya_wcs$site)


## Management column

# Make all entries lower case
catch_kenya_wcs$management <- str_to_lower(catch_kenya_wcs$management)

# Delete management categories with few entries and ambiguous names
catch_kenya_wcs <- catch_kenya_wcs %>%
  filter(!(management %in% c("beachseine", "beachseine and ringnet",
    "ringnet and beach seine", "new community closure")))

# Correct management category names
catch_kenya_wcs$management <- str_replace_all(catch_kenya_wcs$management,
  pattern = "no restriction - mixed gear",
  replacement = "no restriction")

catch_kenya_wcs$management <- str_replace_all(catch_kenya_wcs$management,
  pattern = "pre-community closure",
  replacement = "no restriction")

catch_kenya_wcs$management <- str_replace_all(catch_kenya_wcs$management,
  pattern = "community closure",
  replacement = "closure")

catch_kenya_wcs$management <- str_replace_all(catch_kenya_wcs$management,
  pattern = "old community closure",
  replacement = "closure")

catch_kenya_wcs$management <- str_replace_all(catch_kenya_wcs$management,
  pattern = "old closure",
  replacement = "closure")


## Percent Sampled column

# Fill in 100 where current is na
for(i in 1:nrow(catch_kenya_wcs)){
  
  if(is.na(catch_kenya_wcs$percent_sampled[i]) == TRUE){
    
    catch_kenya_wcs$percent_sampled[i] <- "100"
    
  }
  
}

# Make numeric
catch_kenya_wcs$percent_sampled <- as.numeric(catch_kenya_wcs$percent_sampled)


## Gear column

# Make all entries lower case
catch_kenya_wcs$gear <- str_to_lower(catch_kenya_wcs$gear)

# Correct entries - Basket Trap
catch_kenya_wcs$gear[catch_kenya_wcs$gear == "basket trap"] <- "trap"
catch_kenya_wcs$gear[catch_kenya_wcs$gear == "traps"] <- "trap"

# Correct entries - beachseine
catch_kenya_wcs$gear[catch_kenya_wcs$gear == "beach seine"] <- "beachseine"

# Correct entries - spear
catch_kenya_wcs$gear[catch_kenya_wcs$gear == "diving-spear"] <- "spear"
catch_kenya_wcs$gear[catch_kenya_wcs$gear == "spear diving"] <- "spear"
catch_kenya_wcs$gear[catch_kenya_wcs$gear == "spear gun"] <- "spear"

# Correct entries - gillnet
catch_kenya_wcs$gear[catch_kenya_wcs$gear == "giilnet"] <- "gillnet"
catch_kenya_wcs$gear[catch_kenya_wcs$gear == "gill net"] <- "gillnet"
catch_kenya_wcs$gear[catch_kenya_wcs$gear == "set net"] <- "gillnet"
catch_kenya_wcs$gear[catch_kenya_wcs$gear == "setnet"] <- "gillnet"

# Correct entries - handline
catch_kenya_wcs$gear[catch_kenya_wcs$gear == "hl"] <- "handline"
catch_kenya_wcs$gear[catch_kenya_wcs$gear == "line"] <- "handline"

# Correct entries - monofilament
catch_kenya_wcs$gear[catch_kenya_wcs$gear == "monofillament"] <- "monofilament"

# Correct entries - reefseine
catch_kenya_wcs$gear[catch_kenya_wcs$gear == "reef seine"] <- "reefseine"


## Crew size column

# Make numeric
catch_kenya_wcs$crew_size <- as.numeric(catch_kenya_wcs$crew_size)


## Reshape data

# Make each row contain one individual fish length measurement
catch_kenya_wcs <- pivot_longer(catch_kenya_wcs,
  cols = "1":"32",
  names_to = NULL,
  values_to = "length_cm",
  values_drop_na = TRUE)

# Select only desired columns
catch_kenya_wcs <- select(catch_kenya_wcs, date, year, site, management, crew_id,
  crew_size, gear, percent_sampled, species, length_cm)

# Make length numeric
catch_kenya_wcs$length_cm <- as.numeric(catch_kenya_wcs$length_cm)


## Species column

# Make all scientific names sentence case
catch_kenya_wcs$species <- str_to_sentence(catch_kenya_wcs$species)

# Remove white space
catch_kenya_wcs$species <- str_squish(catch_kenya_wcs$species)

# Find missepllings with taxize package
species_wcs <- unique(catch_kenya_wcs$species)
temp <- gnr_resolve(species_wcs, best_match_only = TRUE, canonical = TRUE)
temp <- filter(temp, temp$user_supplied_name != temp$matched_name2)

# Edit taxize results for Lutjanus lentjan
for (i in 1:nrow(temp)){
  
  if (temp$user_supplied_name[i] == "Lutjanus lentjan"){
    temp$matched_name2[i] <- "Lethrinus lentjan"
  }
  
}

# Replace incorrect/misspelled names
for (i in 1:nrow(temp)){
  catch_kenya_wcs$species <- str_replace_all(catch_kenya_wcs$species,
    temp$user_supplied_name[i], temp$matched_name2[i])
}

# Replace additional incorrect names
catch_kenya_wcs$species <- str_replace_all(catch_kenya_wcs$species,
  pattern = "Synododus variegatus", replacement = "Synodus variegatus")

# Delete unresolvable names (each only occur once)
catch_kenya_wcs <- filter(catch_kenya_wcs, species != "Leptocerus variegatus")
catch_kenya_wcs <- filter(catch_kenya_wcs, species != "Acanthinus trifasciatus")
catch_kenya_wcs <- filter(catch_kenya_wcs, species != "Scavengers")
catch_kenya_wcs <- filter(catch_kenya_wcs, species != "X")
catch_kenya_wcs <- filter(catch_kenya_wcs, species != "Xxx")
catch_kenya_wcs <- filter(catch_kenya_wcs, species != "Xy")
  
# Get WoRMS matches
species_wcs <- unique(catch_kenya_wcs$species)
worms_wcs <- worms_taxa2(species_wcs)

# Correct species names in catch data
for (i in 1:nrow(catch_kenya_wcs)){
  
  # Get index of row in WoRMS results
  x <- which(worms_wcs$supplied_name == catch_kenya_wcs$species[i])
  
  # Replace supplied name with new name
  catch_kenya_wcs$species[i] <- worms_wcs$valid_name[x]
  
}


## Remove rows with unrealistically large fishes

# Remove rows with NA for length
catch_kenya_wcs <- catch_kenya_wcs[!is.na(catch_kenya_wcs$length_cm),]

# Empty vector to hold row numbers to be deleted
temp <- c()

# Loop to add row numbers to temp
for (i in 1:nrow(catch_kenya_wcs)){
  
  # Find index for this species in species_data
  j <- which(species_data$species %in% catch_kenya_wcs$species[i])
  
  # If there was a match
  if(length(j) > 0){
    
    # And if fish length was greater than max length
    if(catch_kenya_wcs$length_cm[i] > species_data$max.length.tl[j]){
      
      # Save row number to temp
      temp <- append(temp, i)
      
    }
    
  }
  
}

# Remove rows with big fishes!
catch_kenya_wcs <- catch_kenya_wcs[-temp, ]


## Add weight column

# Prepare empty column
catch_kenya_wcs$weight_kg <- NA

# Add estimated fish weights
catch_kenya_wcs <- fish_weight(catch_kenya_wcs, key = species_data, worms = worms_wcs)


## Add nutrient weights and L/Lopt

# Nutrient weights
catch_kenya_wcs <- nutrient_weights(catch_kenya_wcs, species_data, worms_wcs)

# L/Lopt
catch_kenya_wcs <- lopt_indicator(catch_kenya_wcs, species_data, worms_wcs)


## Add unique trip id's to each trip

# Change crew size column to effort
catch_kenya_wcs <- rename(catch_kenya_wcs,
  effort = crew_size)

# Unique list of date, site, and crew id combinations
temp <- unique(catch_kenya_wcs[c("date","site","crew_id")])

# Add empty column to hold trip id
temp$trip <- NA

# Add trip id's
for (i in 1:nrow(temp)){
  
  # Add trip id
  temp$trip[i] <- trips[i]
  
}

# Delete all used trip id's from master list
trips <- trips[-(1:i)]

# Add trip id column to catch data
catch_kenya_wcs$trip <- NA

# Add correct trip id's to catch_kenya_wcs
for (i in 1:nrow(catch_kenya_wcs)){
  
  # Filter trip id's to match catch data
  temp2 <- filter(temp, date == catch_kenya_wcs$date[i])
  temp2 <- filter(temp2, site == catch_kenya_wcs$site[i])
  temp2 <- filter(temp2, crew_id == catch_kenya_wcs$crew_id[i])
  
  # Save trip id
  x <- temp2$trip[1]
  catch_kenya_wcs$trip[i] <- x
  
}


## Add missing columns
catch_kenya_wcs$source <- "wcs"
catch_kenya_wcs$country <- "kenya"


## Calculate trip data
trip_kenya_wcs <- fish2trip(catch_kenya_wcs, dri)


## Add missing columns

# Empty columns
trip_kenya_wcs$subsample <- NA

# Add percent sampled to trip data
for (i in 1:nrow(trip_kenya_wcs)){
  
  # Find first occurrence of this trip in catch data
  j <- match(trip_kenya_wcs$trip[i], catch_kenya_wcs$trip)
  
  # Extract percent sampled
  x <- catch_kenya_wcs$percent_sampled[j]
  
  # Save to trip data
  trip_kenya_wcs$subsample[i] <- x
  
}

# Correct weights based on subsamples
trip_kenya_wcs$landings_kg <- trip_kenya_wcs$landings_kg * (100 / trip_kenya_wcs$subsample)
trip_kenya_wcs$calcium_mg <- trip_kenya_wcs$calcium_mg * (100 / trip_kenya_wcs$subsample)
trip_kenya_wcs$iron_mg <- trip_kenya_wcs$iron_mg * (100 / trip_kenya_wcs$subsample)
trip_kenya_wcs$omega.3_g <- trip_kenya_wcs$omega.3_g * (100 / trip_kenya_wcs$subsample)
trip_kenya_wcs$selenium_ug <- trip_kenya_wcs$selenium_ug * (100 / trip_kenya_wcs$subsample)
trip_kenya_wcs$vitamin.a_ug <- trip_kenya_wcs$vitamin.a_ug * (100 / trip_kenya_wcs$subsample)
trip_kenya_wcs$zinc_mg <- trip_kenya_wcs$zinc_mg * (100 / trip_kenya_wcs$subsample)

# Correct yields based on subsamples
trip_kenya_wcs$cpue_kg.effort.day <- trip_kenya_wcs$cpue_kg.effort.day * (100 / trip_kenya_wcs$subsample)
trip_kenya_wcs$calcium_mg.pue <- trip_kenya_wcs$calcium_mg.pue * (100 / trip_kenya_wcs$subsample)
trip_kenya_wcs$iron_mg.pue <- trip_kenya_wcs$iron_mg.pue * (100 / trip_kenya_wcs$subsample)
trip_kenya_wcs$omega.3_g.pue <- trip_kenya_wcs$omega.3_g.pue * (100 / trip_kenya_wcs$subsample)
trip_kenya_wcs$selenium_ug.pue <- trip_kenya_wcs$selenium_ug.pue * (100 / trip_kenya_wcs$subsample)
trip_kenya_wcs$vitamin.a_ug.pue <- trip_kenya_wcs$vitamin.a_ug.pue * (100 / trip_kenya_wcs$subsample)
trip_kenya_wcs$zinc_mg.pue <- trip_kenya_wcs$zinc_mg.pue * (100 / trip_kenya_wcs$subsample)
trip_kenya_wcs$nutrients_pdv.pue <- ((trip_kenya_wcs$nutrients_pdv * 10) * trip_kenya_wcs$cpue_kg.effort.day) / 6

# Correct pdv yields based on corrected yields
trip_kenya_wcs$calcium_pdv.pue <- trip_kenya_wcs$calcium_mg.pue / dri[dri$nutrients == "calcium_mg", "dri"]
trip_kenya_wcs$iron_pdv.pue <- trip_kenya_wcs$iron_mg.pue / dri[dri$nutrients == "iron_mg", "dri"]
trip_kenya_wcs$omega.3_pdv.pue <- trip_kenya_wcs$omega.3_g.pue / dri[dri$nutrients == "omega3_g", "dri"]
trip_kenya_wcs$selenium_pdv.pue <- trip_kenya_wcs$selenium_ug.pue / dri[dri$nutrients == "selenium_ug", "dri"]
trip_kenya_wcs$vitamin.a_pdv.pue <- trip_kenya_wcs$vitamin.a_ug.per.100g / dri[dri$nutrients == "vitamina_ug", "dri"]
trip_kenya_wcs$zinc_pdv.pue <- trip_kenya_wcs$zinc_mg.pue / dri[dri$nutrients == "zinc_mg", "dri"]

# Make date column as date
trip_kenya_wcs$date <- as.Date(trip_kenya_wcs$date, origin = "1970-01-01")


## Add to catch_kenya

# Aggregate by landing site and species
temp <- aggregate(weight_kg ~ species + site,
  data = catch_kenya_wcs, FUN = sum)

# Make a summary of landings by species by site
temp <- temp %>%
  pivot_wider(names_from = site, values_from = weight_kg)

# Bind to catch_kenya master
catch_kenya <- full_join(catch_kenya, temp)


## Add to trip data
trip <- bind_rows(trip, trip_kenya_wcs)
trip <- trip %>% select(-effort, -subsample)


## Add to fish data
fish <- bind_rows(fish, catch_kenya_wcs)
fish <- fish %>% select(-(year:zinc_mg))




##### Blue Ventures #####

## This section cleans data provided by Blue Ventures in Madagascar

# NB: Ampandikoara and following are all in grams and need to be converted
# to kg. At Ankify, it turns back to kg.

# Load data
catch_mdg_bv <- read_csv("data/raw-data/catch-monitoring/madagascar-catch-monitoring-blueventures.csv", 
    col_types = cols(`2010` = col_number(), 
        `2011` = col_number(), `2012` = col_number(), 
        `2013` = col_number(), `2014` = col_number(), 
        `2015` = col_number(), `2016` = col_number(), 
        `2017` = col_number(), `2018` = col_number(), 
        `2019` = col_number(), `2020` = col_number(), 
        `2021` = col_number(), `2022` = col_number()))

# Correct column names
catch_mdg_bv <- catch_mdg_bv %>% rename(c(
  site = "Landing_site",
  species = "Species"
))

# Make landing site names lower case
catch_mdg_bv$site <- str_to_lower(catch_mdg_bv$site)


## Correct species names

# Remove all potentially confounding characters from species list
catch_mdg_bv$species <- str_replace_all(catch_mdg_bv$species,
  pattern = "\\(",
  replacement = "")
catch_mdg_bv$species <- str_replace_all(catch_mdg_bv$species,
  pattern = "\\)",
  replacement = "")
catch_mdg_bv$species <- str_replace_all(catch_mdg_bv$species,
  pattern = "<",
  replacement = "")
catch_mdg_bv$species <- str_replace_all(catch_mdg_bv$species,
  pattern = ">",
  replacement = "")
catch_mdg_bv$species <- str_replace_all(catch_mdg_bv$species,
  pattern = "\\.",
  replacement = "")
catch_mdg_bv$species <- str_replace_all(catch_mdg_bv$species,
  pattern = "�",
  replacement = " ")

# Manually replace known problems
catch_mdg_bv$species <- str_replace_all(catch_mdg_bv$species,
  pattern = "Fistularidae",
  replacement = "Fistulariidae")
catch_mdg_bv$species <- str_replace_all(catch_mdg_bv$species,
  pattern = "Blennidae",
  replacement = "Blenniidae")
catch_mdg_bv$species <- str_replace_all(catch_mdg_bv$species,
  pattern = "Carcharinidae",
  replacement = "Carcharhinidae")
catch_mdg_bv$species <- str_replace_all(catch_mdg_bv$species,
  pattern = "Clupidae",
  replacement = "Clupeidae")
catch_mdg_bv$species <- str_replace_all(catch_mdg_bv$species,
  pattern = "Hemirhamphidae",
  replacement = "Hemiramphidae")
catch_mdg_bv$species <- str_replace_all(catch_mdg_bv$species,
  pattern = "Muglidae",
  replacement = "Mugilidae")
catch_mdg_bv$species <- str_replace_all(catch_mdg_bv$species,
  pattern = "Mulloidicthys",
  replacement = "Mulloidichthys")
catch_mdg_bv$species <- str_replace_all(catch_mdg_bv$species,
  pattern = "Taenuria melanospilos",
  replacement = "Taeniura melanospilos")
catch_mdg_bv$species <- str_replace_all(catch_mdg_bv$species,
  pattern = "Charcharinus sp",
  replacement = "Carcharhinus")
catch_mdg_bv$species <- str_replace_all(catch_mdg_bv$species,
  pattern = "Ppalunirus sp",
  replacement = "Palinurus")
catch_mdg_bv$species <- str_replace_all(catch_mdg_bv$species,
  pattern = "Carcharinus",
  replacement = "Carcharhinus")
catch_mdg_bv$species <- str_replace_all(catch_mdg_bv$species,
  pattern = "Rhinbatus sp",
  replacement = "Rhinobatos")
catch_mdg_bv$species <- str_replace_all(catch_mdg_bv$species,
  pattern = "Caesion sp",
  replacement = "Caesio")
catch_mdg_bv$species <- str_replace_all(catch_mdg_bv$species,
  pattern = "Carcharininus sp",
  replacement = "Carcharhinus")
catch_mdg_bv$species <- str_replace_all(catch_mdg_bv$species,
  pattern = "Sphyrnea sp",
  replacement = "Sphyraena")

# Remove unidentified species
catch_mdg_bv <- filter(catch_mdg_bv, species != "TBI")
catch_mdg_bv <- filter(catch_mdg_bv, species != "mixed species")
catch_mdg_bv <- filter(catch_mdg_bv, species != "mixed species 5cm")

# Find misspellings with taxize package
species_bv <- unique(catch_mdg_bv$species)
temp <- gnr_resolve(species_bv, best_match_only = TRUE, canonical = TRUE)
temp <- filter(temp, temp$user_supplied_name != temp$matched_name2)

# Replace misspellings
for (i in 1:nrow(temp)){
  
  catch_mdg_bv$species <- str_replace_all(catch_mdg_bv$species,
    pattern = temp$user_supplied_name[i],
    replacement = temp$matched_name2[i])
  
}

# Get WoRMS matches
species_bv <- unique(catch_mdg_bv$species)
worms_bv <- worms_taxa2(species_bv)

# Confirm all species were matched by worms
temp <- species_bv[!(species_bv %in% worms_bv$supplied_name)]

# Correct species names in catch data
for (i in 1:nrow(catch_mdg_bv)){
  
  # Get index of row in WoRMS results
  x <- which(worms_bv$supplied_name == catch_mdg_bv$species[i])
  
  # Replace supplied name with new name
  catch_mdg_bv$species[i] <- worms_bv$valid_name[x]
  
}

# Remove non-chordata
temp <- subset(worms_bv$valid_name, worms_bv$phylum != "Chordata")
catch_mdg_bv <- filter(catch_mdg_bv, !(catch_mdg_bv$species %in% temp))
worms_bv <- filter(worms_bv, phylum == "Chordata")


## Convert landing sites in g to kg

# Vector of all sites
sites <- unique(catch_mdg_bv$site)

# Vector of sites in grams
sites2 <- sites[which(sites == "ampandikoara") : which(sites == "lovobe")]

# Convert to kilograms
catch_mdg_bv[catch_mdg_bv$site %in% sites2, 3:15] <- catch_mdg_bv[catch_mdg_bv$site %in% sites2, 3:15] / 1000


## Add to master data frame

# Aggregate by species and site

# Aggregate each year individually
temp01 <- aggregate(`2010` ~ site + species,
  data = catch_mdg_bv, FUN = sum, na.rm = TRUE)
temp02 <- aggregate(`2011` ~ site + species,
  data = catch_mdg_bv, FUN = sum, na.rm = TRUE)
temp03 <- aggregate(`2012` ~ site + species,
  data = catch_mdg_bv, FUN = sum, na.rm = TRUE)
temp04 <- aggregate(`2013` ~ site + species,
  data = catch_mdg_bv, FUN = sum, na.rm = TRUE)
temp05 <- aggregate(`2014` ~ site + species,
  data = catch_mdg_bv, FUN = sum, na.rm = TRUE)
temp06 <- aggregate(`2015` ~ site + species,
  data = catch_mdg_bv, FUN = sum, na.rm = TRUE)
temp07 <- aggregate(`2016` ~ site + species,
  data = catch_mdg_bv, FUN = sum, na.rm = TRUE)
temp08 <- aggregate(`2017` ~ site + species,
  data = catch_mdg_bv, FUN = sum, na.rm = TRUE)
temp09 <- aggregate(`2018` ~ site + species,
  data = catch_mdg_bv, FUN = sum, na.rm = TRUE)
temp10 <- aggregate(`2019` ~ site + species,
  data = catch_mdg_bv, FUN = sum, na.rm = TRUE)
temp11 <- aggregate(`2020` ~ site + species,
  data = catch_mdg_bv, FUN = sum, na.rm = TRUE)
temp12 <- aggregate(`2021` ~ site + species,
  data = catch_mdg_bv, FUN = sum, na.rm = TRUE)
temp13 <- aggregate(`2022` ~ site + species,
  data = catch_mdg_bv, FUN = sum, na.rm = TRUE)

# Combine into dataframe
temp <- full_join(temp01, temp02, by = c("species", "site"))
temp <- full_join(temp, temp03, by = c("species", "site"))
temp <- full_join(temp, temp04, by = c("species", "site"))
temp <- full_join(temp, temp05, by = c("species", "site"))
temp <- full_join(temp, temp06, by = c("species", "site"))
temp <- full_join(temp, temp07, by = c("species", "site"))
temp <- full_join(temp, temp08, by = c("species", "site"))
temp <- full_join(temp, temp09, by = c("species", "site"))
temp <- full_join(temp, temp10, by = c("species", "site"))
temp <- full_join(temp, temp11, by = c("species", "site"))
temp <- full_join(temp, temp12, by = c("species", "site"))
temp <- full_join(temp, temp13, by = c("species", "site"))
catch_mdg_bv <- temp

# Add a column for row means
temp <- catch_mdg_bv %>% select("2010":"2022")
temp <- data.matrix(temp)
catch_mdg_bv$mean <- rowMeans(temp, na.rm = TRUE)

# Remove year columns
catch_mdg_bv <- catch_mdg_bv %>% select(species, site, mean)

# Aggregate by landing site and species to remove duplicates
catch_mdg_bv <- aggregate(mean ~ species + site,
  data = catch_mdg_bv, FUN = sum)

# Make a summary of landings by species by site
catch_mdg_bv <- catch_mdg_bv %>%
  pivot_wider(names_from = site, values_from = mean)

# Replace NA with 0
catch_mdg_bv[is.na(catch_mdg_bv)] <- 0

# Save catch_mdg
catch_mdg <- catch_mdg_bv




##### Tiliouine 2019 #####

## This section cleans and compiles data from:

## Tiliouine, O. (2019). An investigation into fisheries trends and a critical analysis
##    of fisheries management within the Bay of Ranobe, Madagascar. University College 
##    London.

## NB: The observations in this data are aggregated by species, so each row includes
## weight and length measurements for one fish as well as the total weight of all
## individuals of that species in the catch. I am importing the length and total weight,
## but not the subsampled weight. This means that length and weight will not reflect
## known weight-length relationships.

# Import data
df1 <- read_excel("data/raw-data/catch-monitoring/madagascar-catch-monitoring-tiliouine2019a.xlsx")
df2 <- read_excel("data/raw-data/catch-monitoring/madagascar-catch-monitoring-tiliouine2019b.xlsx")
df3 <- read_excel("data/raw-data/catch-monitoring/madagascar-catch-monitoring-tiliouine2019c.xlsx")
df4 <- read_excel("data/raw-data/catch-monitoring/madagascar-catch-monitoring-tiliouine2019d.xlsx")

# Select and rename columns
df1 <- df1 %>% select("Village...2", "Date...3", "Ticket...4", "Species", "Tw (kg)",
  "TL (cm)", "Gear 1", "Nb pers")
df1 <- df1 %>% rename(site = "Village...2",
  date = "Date...3",
  crew_id = "Ticket...4",
  species = "Species",
  weight_kg = "Tw (kg)",
  length_cm = "TL (cm)",
  gear = "Gear 1",
  effort = "Nb pers")

df2 <- df2 %>% select("Village...2", "Date...3", "Ticket...4", "Species", "Tw(kg)",
  "L (cm)", "Gear", "nb.pers")
df2 <- df2 %>% rename(site = "Village...2",
  date = "Date...3",
  crew_id = "Ticket...4",
  species = "Species",
  weight_kg = "Tw(kg)",
  length_cm = "L (cm)",
  gear = "Gear",
  effort = "nb.pers")

df3 <- df3 %>% select("Village...2", "Date...3", "Ticket...4", "Species", "Tw(kg)",
  "L (cm)", "Gear", "nb.pers")
df3 <- df3 %>% rename(site = "Village...2",
  date = "Date...3",
  crew_id = "Ticket...4",
  species = "Species",
  weight_kg = "Tw(kg)",
  length_cm = "L (cm)",
  gear = "Gear",
  effort = "nb.pers")

df4 <- df4 %>% select("Village...2", "Date...3", "Ticket...4", "Species",
  "Total weight (kg)", "TL (cm)", "Gear 1", "nb.pers")
df4 <- df4 %>% rename(site = "Village...2",
  date = "Date...3",
  crew_id = "Ticket...4",
  species = Species,
  weight_kg = "Total weight (kg)",
  length_cm = "TL (cm)",
  gear = "Gear 1",
  effort = "nb.pers")

# Set column types
df1$site <- as.character(df1$site)
df1$date <- as.character(df1$date)
df1$crew_id <- as.numeric(df1$crew_id)
df1$species <- as.character(df1$species)
df1$weight_kg <- as.numeric(df1$weight_kg)
df1$length_cm <- as.numeric(df1$length_cm)
df1$gear <- as.character(df1$gear)
df1$effort <- as.numeric(df1$effort)

df2$site <- as.character(df2$site)
df2$date <- as.character(df2$date)
df2$crew_id <- as.numeric(df2$crew_id)
df2$species <- as.character(df2$species)
df2$weight_kg <- as.numeric(df2$weight_kg)
df2$length_cm <- as.numeric(df2$length_cm)
df2$gear <- as.character(df2$gear)
df2$effort <- as.numeric(df2$effort)

df3$site <- as.character(df3$site)
df3$date <- as.character(df3$date)
df3$crew_id <- as.numeric(df3$crew_id)
df3$species <- as.character(df3$species)
df3$weight_kg <- as.numeric(df3$weight_kg)
df3$length_cm <- as.numeric(df3$length_cm)
df3$gear <- as.character(df3$gear)
df3$effort <- as.numeric(df3$effort)

df4$site <- as.character(df4$site)
df4$date <- as.character(df4$date)
df4$crew_id <- as.numeric(df4$crew_id)
df4$species <- as.character(df4$species)
df4$weight_kg <- as.numeric(df4$weight_kg)
df4$length_cm <- as.numeric(df4$length_cm)
df4$gear <- as.character(df4$gear)
df4$effort <- as.numeric(df4$effort)

# Combine years
catch_mdg_tiliouine <- bind_rows(df1, df2, df3, df4)


## Site column

# Make lower case
catch_mdg_tiliouine$site <- str_to_lower(catch_mdg_tiliouine$site)

# Ambolomailaky
catch_mdg_tiliouine$site <- catch_mdg_tiliouine$site %>%
  str_replace_all(pattern = "ambolimailaka", replacement = "temp")
catch_mdg_tiliouine$site <- catch_mdg_tiliouine$site %>%
  str_replace_all(pattern = "ambolimailaky", replacement = "temp")
catch_mdg_tiliouine$site <- catch_mdg_tiliouine$site %>%
  str_replace_all(pattern = "ambolomailaky", replacement = "temp")
catch_mdg_tiliouine$site <- catch_mdg_tiliouine$site %>%
  str_replace_all(pattern = "amb", replacement = "temp")
catch_mdg_tiliouine$site <- catch_mdg_tiliouine$site %>%
  str_replace_all(pattern = "temp", replacement = "ambolomailaky")

# Andrevo
catch_mdg_tiliouine$site <- catch_mdg_tiliouine$site %>%
  str_replace_all(pattern = "andrevo", replacement = "temp")
catch_mdg_tiliouine$site <- catch_mdg_tiliouine$site %>%
  str_replace_all(pattern = "and", replacement = "temp")
catch_mdg_tiliouine$site <- catch_mdg_tiliouine$site %>%
  str_replace_all(pattern = "an", replacement = "temp")
catch_mdg_tiliouine$site <- catch_mdg_tiliouine$site %>%
  str_replace_all(pattern = "temp", replacement = "andrevo")

# Beravy
catch_mdg_tiliouine$site <- catch_mdg_tiliouine$site %>%
  str_replace_all(pattern = "baravy", replacement = "temp")
catch_mdg_tiliouine$site <- catch_mdg_tiliouine$site %>%
  str_replace_all(pattern = "beravy", replacement = "temp")
catch_mdg_tiliouine$site <- catch_mdg_tiliouine$site %>%
  str_replace_all(pattern = "be", replacement = "temp")
catch_mdg_tiliouine$site <- catch_mdg_tiliouine$site %>%
  str_replace_all(pattern = "bl", replacement = "temp")
catch_mdg_tiliouine$site <- catch_mdg_tiliouine$site %>%
  str_replace_all(pattern = "temp", replacement = "beravy")

# Ifaty
catch_mdg_tiliouine$site <- catch_mdg_tiliouine$site %>%
  str_replace_all(pattern = "ifaty", replacement = "temp")
catch_mdg_tiliouine$site <- catch_mdg_tiliouine$site %>%
  str_replace_all(pattern = "if", replacement = "temp")
catch_mdg_tiliouine$site <- catch_mdg_tiliouine$site %>%
  str_replace_all(pattern = "temp", replacement = "ifaty")

# Delete empty rows
catch_mdg_tiliouine <- catch_mdg_tiliouine[-c(2494, 5704, 8840, 10201),]


## Date column

# Make a character vector
catch_mdg_tiliouine$date <- as.character(catch_mdg_tiliouine$date)

# The first format to fix - dd.mm.yyyy
v1 <- catch_mdg_tiliouine$date[str_which(catch_mdg_tiliouine$date, pattern = "\\.")]

# Empty date vector
v2 <- vector(mode = "character", length = length(v1))

# List object for all of these dates
temp <- str_split(v1, "\\.")

# Re-assemble dates to a better format
for (i in 1:length(v1)){
  
  v2[i] <- paste("2014", temp[i][[1]][2], temp[i][[1]][1], sep = "-")
  
}

# Save as date
v2 <- as.Date(v2)

# The second format to fix - "#####"
v3 <- catch_mdg_tiliouine$date[str_length(catch_mdg_tiliouine$date) == 5]

# Make numeric
v3 <- as.numeric(v3)

# Make date
v3 <- as.Date(v3, origin = "1899-12-31")

# Extract all other dates
v4 <- catch_mdg_tiliouine$date[str_which(catch_mdg_tiliouine$date, pattern = "-")]

# Make date
v4 <- as.Date(v4)

# Make one vector
v5 <- c(v2, v3, v4)

# Replace date column
catch_mdg_tiliouine$date <- v5


## Species column

# List of unique species
species_tiliouine <- unique(catch_mdg_tiliouine$species)

# List of matching scientific names (for those that can be matched)
scinames <- comm2sci(species_tiliouine, db = "worms")

# List of unambiguous scientific name matches
scinames2 <- scinames[lengths(scinames) == 1]

# Turn list into dataframe
scinames2 <- bind_rows(scinames2)
scinames2 <- pivot_longer(scinames2, everything(), names_to = "common", values_to = "scientific")

# Replace applicable names in catch data
for (i in 1:nrow(scinames2)){
  
  catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == scinames2$common[i]] <-
    scinames2$scientific[i]
  
}

# A list of unmatched names
comnames <- names(scinames[lengths(scinames) != 1])

# Dataframe of possible matches from fishbase
scinames <- common_to_sci(comnames)

# Dataframe of species with only one match (unambiguous matches)
scinames2 <- scinames %>%
  group_by(ComName) %>%
  filter(n() == 1)

# Replace unambiguous matched names in catch data
for(i in 1:nrow(scinames2)){
  
  catch_mdg_tiliouine$species[
    grabl(catch_mdg_tiliouine$species, pattern = scinames2$ComName[i],
      maxDist = 2)
  ] <- scinames2$Species[i]
  
}

# Manually replace remaining common names
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Whitespotted Pufferfish"] <- "Arothron meleagris"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Dark Damsel fish"] <- "Pycnochromis leucurus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "White-spotted guitarfish"] <- "Rhynchobatus djiddensis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Afr. W. Spot rabbitfish"] <- "Siganus sutor"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Redgill emperor"] <- "Lethrinus rubrioperculatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Ring-tailed cardinalfish"] <- "Ostorhinchus aureus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Dark damselfish"] <- "Pycnochromis leucurus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "A W Rabbit fish"] <- "Siganus sutor"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "A W Rabbitfish"] <- "Siganus sutor"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "A. W. Rabbit Fish"] <- "Siganus sutor"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "A.W Rabbitfish"] <- "Siganus sutor"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "A.W. Rabbit Fish"] <- "Siganus sutor"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "A.W. Rabbitfish"] <- "Siganus sutor"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "a.w. spot rabbitfish"] <- "Siganus sutor"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Abalily (Mix)"] <- "other"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Afr. W. spot rabbitfish"] <- "Siganus sutor"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Afr.W. Spot Rabbitfish"] <- "Siganus sutor"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Afr.W.Spot Rabbitfish"] <- "Siganus sutor"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "African white spotted rabbit fish"] <- "Siganus sutor"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "African white spotted rabbitfish"] <- "Siganus sutor"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "african whitespotted rabbitfish"] <- "Siganus sutor"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "African whitespotted rabbitfish"] <- "Siganus sutor"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "African Whitespotted rabbitfish"] <- "Siganus sutor"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "African whrite spotted rabbit fish"] <- "Siganus sutor"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "African whrite stoffet rabbit"] <- "Siganus sutor"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Ambalily"] <- "other"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "AmBalily"] <- "other"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Ambalily (mix)"] <- "other"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Ambotsoky"] <- "Sillago sihama"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Anchory spp"] <- "Engraulidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Anchovy spp."] <- "Engraulidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Anchovy Spp."] <- "Engraulidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Arabian spine cheek"] <- "Scolopsis ghanam"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Arabian spinecheek"] <- "Scolopsis ghanam"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Arabian Spinecheek"] <- "Scolopsis ghanam"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Arabianspine cheeck"] <- "Scolopsis ghanam"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Barfin moray"] <- "Gymnothorax zonipectis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Barracuda spp."] <- "Sphyraena"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Bartail flat herd"] <- "Platycephalus indicus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Bat fish"] <- "Ogcocephalidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Batfish"] <- "Ogcocephalidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Batfish spp."] <- "Ogcocephalidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Beauferts crocodilefish"] <- "other"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Betampy"] <- "other"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Bicolor parrot fish"] <- "Scaridae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Big eye emperor"] <- "Monotaxis grandoculis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Big eye trevally"] <- "Caranx sexfasciatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Big Eye trevally"] <- "Caranx sexfasciatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Big-eye Emperor"] <- "Monotaxis grandoculis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Bigeye emperor"] <- "Monotaxis grandoculis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Bigeye Trevally spp."] <- "Caranx sexfasciatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Bigeye Emperor"] <- "Monotaxis grandoculis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Black barred parrotfish"] <- "Scaridae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Black fin squirrel fish"] <- "Neoniphon opercularis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "black nape large eye bream"] <- "Gymnocranius"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Black nape large eye bream"] <- "Gymnocranius"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Black Peak Surgeon Fish"] <- "Acanthurus nigricauda"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Black Spot Emperor"] <- "Lethrinus harak"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Black spot snapper"] <- "Lutjanus monostigma"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Black Spot Snapper"] <- "Lutjanus monostigma"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Black spotted butter"] <- "other"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Black spotted electric ray"] <- "Torpedo fuscomaculata"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Black Spotted Sweetlips"] <- "Plectorhinchus gibbosus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Black streak sureonfish"] <- "Acanthurus nigricauda"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Black streak surgeon fish"] <- "Acanthurus nigricauda"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Black Streak Surgeonfish"] <- "Acanthurus nigricauda"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Black stripe cardinal fish"] <- "Ostorhinchus nigrofasciatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Black striped goat fish"] <- "Upeneus tragula"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Black stwent surseen fish"] <- "Acanthurus nigricauda"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Black-spot snapper"] <- "Lutjanus monostigma"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Black-spot Snapper"] <- "Lutjanus monostigma"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Black-spotted sweetlips"] <- "Plectorhinchus gibbosus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Blackbanded cardinalfish"] <- "Ostorhinchus nigrofasciatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Blackfin porcupinefish"] <- "Diodontidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "blacknap large eye bream"] <- "Gymnocranius"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Blackspot emporer"] <- "Lethrinus harak"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Blackspot snapper"] <- "Lutjanus monostigma"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Blackspot Snapper"] <- "Lutjanus monostigma"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Blackspotted electricray"] <- "Torpedo fuscomaculata"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "blackspotted sweetlips"] <- "Plectorhinchus gibbosus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Blackspotted sweetlips"] <- "Plectorhinchus gibbosus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Blackspotted Sweetlips"] <- "Plectorhinchus gibbosus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Blackstreak surgeon"] <- "Acanthurus nigricauda"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Blackstreak surgeonfish"] <- "Acanthurus nigricauda"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Blackstreak Surgeonfish"] <- "Acanthurus nigricauda"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Black streak surgeonfish"] <- "Acanthurus nigricauda"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Blod Spot Squirrelfish"] <- "Neoniphon sammara"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Blood Spot Squirrelfish"] <- "Neoniphon sammara"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Bloodspot Squirrelfish"] <- "Neoniphon sammara"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Bludger Travelly"] <- "Carangoides gymnostethus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Blue barred paroi fish"] <- "Scarus ghobban"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Blue barred parroi fish"] <- "Scarus ghobban"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Blue barred parroit fish"] <- "Scarus ghobban"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "blue barred parrot fish"] <- "Scarus ghobban"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Blue barred parrot fish"] <- "Scarus ghobban"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "blue barred parrotfish"] <- "Scarus ghobban"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Blue barred parrotfish"] <- "Scarus ghobban"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Blue fin trevally"] <- "Caranx melampygus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Blue spine unicorn"] <- "Naso unicornis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Blue spine unicorn fish"] <- "Naso unicornis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Bluebarred parrotfish"] <- "Scarus ghobban"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Bluebarred Parrotfish"] <- "Scarus ghobban"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Bluefin Trevally spp."] <- "Caranx melampygus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Bluefin Trevelly"] <- "Caranx melampygus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "bluespined unicornfish"] <- "Naso unicornis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Bluespotted stingray"] <- "Neotrygon kuhlii"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Boat fish"] <- "Ogcocephalidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Bohadschia spp."] <- "Bohadschia"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Bonito"] <- "Scombridae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Borea (mix)"] <- "other"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Borea (Mix)"] <- "other"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Borea mix"] <- "other"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Borea Mix"] <- "other"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Broad head flat head"] <- "Sunagocia arenicola"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Broad hehd flutehead"] <- "Sunagocia arenicola"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Bronze soldierfish"] <- "Myripristis adusta"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Bronze Soldierfish"] <- "Myripristis adusta"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Brushtail tang"] <- "Zebrasoma scopas"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Brushtail Tang"] <- "Zebrasoma scopas"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Brustail Tang"] <- "Zebrasoma scopas"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Buffet"] <- "other"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Butter fly fish"] <- "Chaetodontidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Cardinal fish spp."] <- "Apogonidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Cardinalfish spp."] <- "Apogonidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Cathefish"] <- "other"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Catle fish"] <- "other"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Catlefish"] <- "other"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Cheek Board Wrasse"] <- "Halichoeres hortulanus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Cigar wrass"] <- "Cheilio inermis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Circular batfish"] <- "Platax orbicularis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Common Gunner Helmet"] <- "Dactyloptena orientalis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Common Helmet gunard"] <- "Dactyloptena orientalis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Common helmet gurnard"] <- "Dactyloptena orientalis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Common Helmet gurnard"] <- "Dactyloptena orientalis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Common Helmet Gurnard"] <- "Dactyloptena orientalis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Common hilmet guarnard"] <- "Dactyloptena orientalis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "CorGill net fish"] <- "other"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Coris gaimard africana"] <- "Coris gaimard"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Cornet fish"] <- "Fistularia commersonii"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Cornet Fish"] <- "Fistularia commersonii"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Crescent Banded Grunter"] <- "Terapon jarbua"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Crescent wrasse"] <- "Thalassoma lunare"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Cresent wrasse"] <- "Thalassoma lunare"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Cresent-Banded Grunter"] <- "Terapon jarbua"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Creslent wrasse"] <- "Thalassoma lunare"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Crown Squirrel Fish"] <- "Sargocentron diadema"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Dandezo"] <- "other"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Dark damsel"] <- "Pomacentrus aquilus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Dark Damsel"] <- "Pomacentrus aquilus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Dark damsel fish"] <- "Pomacentrus aquilus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Dark Damsel Fish"] <- "Pomacentrus aquilus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Dark dansel"] <- "Pomacentrus aquilus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Dash and Dot Fish"] <- "Parupeneus barberinus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Dash and dot gait fish"] <- "Parupeneus barberinus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Dash and Dot Goat Fish"] <- "Parupeneus barberinus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Dash-and-dot goat fish"] <- "Parupeneus barberinus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Desjarde sailfin tang"] <- "Zebrasoma desjardinii"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Desjarden's sailfin tang"] <- "Zebrasoma desjardinii"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "desjardin sailfin tang"] <- "Zebrasoma desjardinii"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "DesJardin Sailfin Tang"] <- "Zebrasoma desjardinii"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Desjardin's sailfin tang"] <- "Zebrasoma desjardinii"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Desjardin’s sailfin tang"] <- "Zebrasoma desjardinii"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Dragon Wrasse"] <- "Novaculichthys taeniourus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Dusky parroit fish"] <- "Scarus niger"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Dusky surgeon fish"] <- "Acanthurus nigrofuscus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "dusky surgeonfish"] <- "Acanthurus nigrofuscus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Dusky surgeonfish"] <- "Acanthurus nigrofuscus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Dusky Surgeonfish"] <- "Acanthurus nigrofuscus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Dusky sweet lips"] <- "Plectorhinchus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Dusky sweetlips"] <- "Plectorhinchus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Dusky Sweetlips"] <- "Plectorhinchus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "elegant surgeonfish"] <- "Naso elegans"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Emperor spp."] <- "Lethrinidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "False eye sergeant fish"] <- "Abudefduf sparoides"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "False eye sergeantfish"] <- "Abudefduf sparoides"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "False eye sergent"] <- "Abudefduf sparoides"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "File fish spp."] <- "Monacanthidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Filefish"] <- "Monacanthidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Filefish spp."] <- "Monacanthidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Five  lined cardinalfish"] <- "Cheilodipterus quinquelineatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Five line cardinal fish"] <- "Cheilodipterus quinquelineatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Five linecardinal fish"] <- "Cheilodipterus quinquelineatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Five lined cardinal fish"] <- "Cheilodipterus quinquelineatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Flathead spp."] <- "Platycephalidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Flounder spp."] <- "Poecilopsettidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Flying fish"] <- "Exocoetidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Flying Fish spp."] <- "Exocoetidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Flyingfish spp."] <- "Exocoetidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Fork tail rabbit fish"] <- "Siganus argenteus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Geba Bato"] <- "Herklotsichthys quadrimaculatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Geba mena"] <- "Amblygaster sirm"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Geba mena or bato"] <- "Clupeidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Geoelectric moray"] <- "Muraenidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Giant Moray eel"] <- "Gymnothorax javanicus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Giant Trevally spp."] <- "Caranx ignobilis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Giant wuarly"] <- "Osphronemus goramy"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Glass eye big eye"] <- "Priacanthus blochii"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Glass-eye Big Eye"] <- "Priacanthus blochii"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Glass-eye Big-eye"] <- "Priacanthus blochii"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "glasseye bigeye"] <- "Priacanthus blochii"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Glasseye bigeye"] <- "Priacanthus blochii"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Glasseye Bigeye"] <- "Priacanthus blochii"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Glasseye bigeyefish"] <- "Priacanthus blochii"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Goggle-eye bigeye"] <- "Priacanthus blochii"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Gogle-eye bigeye"] <- "Priacanthus blochii"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Gold Barre Wrasse"] <- "Thalassoma hebraicum"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Gold Spot Herring"] <- "Herklotsichthys quadrimaculatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Gold spotted sweet lips"] <- "Plectorhinchus flavomaculatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Gold Spotted Sweet Lips"] <- "Plectorhinchus flavomaculatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Gold spotted sweetlips"] <- "Plectorhinchus flavomaculatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Gold Spotted Sweetlips"] <- "Plectorhinchus flavomaculatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Gold-spotted sweetlips"] <- "Plectorhinchus flavomaculatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Gold-spotted Sweetlips"] <- "Plectorhinchus flavomaculatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Goldlined seabream"] <- "Rhabdosargus sarba"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Goldlined seabream"] <- "Rhabdosargus sarba"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Goldlined seabream"] <- "Rhabdosargus sarba"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Goldring bristletooth"] <- "Ctenochaetus cyanocheilus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "goldspot sweetlips"] <- "Plectorhinchus flavomaculatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Goldspotted sweetlips"] <- "Plectorhinchus flavomaculatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Goldspotted Sweetlips"] <- "Plectorhinchus flavomaculatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Golospotted sweetlips"] <- "Plectorhinchus flavomaculatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Google eye big eye"] <- "Priacanthus blochii"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Google Eye Big Eye"] <- "Priacanthus blochii"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Google-eye Big Eye"] <- "Priacanthus blochii"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Graceful lizand fish"] <- "Saurida gracilis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Graceful lizard fish"] <- "Saurida gracilis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Graceful Lizard Fish"] <- "Saurida gracilis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "graceful lizardfish"] <- "Saurida gracilis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Graceful lizardfish"] <- "Saurida gracilis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Graceful Lizardfish"] <- "Saurida gracilis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Gracefull Lizardfish"] <- "Saurida gracilis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Grouper"] <- "Serranidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Grouper spp."] <- "Serranidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Grunge-lined cardinafish"] <- "Taeniamia fucata"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Hatokantendro"] <- "Scomberoides commersonianus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Herring Spp (Mix)"] <- "Clupeidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Herring spp."] <- "Clupeidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Herring spp. (Mix)"] <- "Clupeidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "High fin rudderfish"] <- "Kyphosus cinerascens"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Highfin rudder fish"] <- "Kyphosus cinerascens"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Highfin rudderfish"] <- "Kyphosus cinerascens"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Highfin Rudderfish"] <- "Kyphosus cinerascens"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Honey comb moray"] <- "Gymnothorax favagineus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Honeycomb sting ray"] <- "Himantura uarnak"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Horse face unicorn fish"] <- "Naso"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Humb back snapper"] <- "Lutjanus gibbus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Humb back unicorn fish"] <- "Naso brachycentron"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Hump back snapper"] <- "Lutjanus gibbus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "I.D. Crocoile fish"] <- "Papilloculiceps longiceps"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "I.O Crocodile fish"] <- "Papilloculiceps longiceps"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "I.O. crocodilefish"] <- "Papilloculiceps longiceps"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "I.O. Crocodilefish"] <- "Papilloculiceps longiceps"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "I.O.crocodilefish"] <- "Papilloculiceps longiceps"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Icillar Half Beak"] <- "Hyporhamphus affinis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Idian ocean crocodilefish"] <- "Papilloculiceps longiceps"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Indian goat fish"] <- "Parupeneus indicus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Indian ocean crocodile fish"] <- "Papilloculiceps longiceps"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Indian Ocean crocodile fish"] <- "Papilloculiceps longiceps"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Indian ocean crocodilefish"] <- "Papilloculiceps longiceps"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Indian Ocean crocodilefish"] <- "Papilloculiceps longiceps"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Indian ocean goat fish"] <- "Parupeneus indicus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Indian oceanne crocodile fish"] <- "Papilloculiceps longiceps"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Indian thead fin"] <- "Leptomelanosoma indicum"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Indian threadfin"] <- "Leptomelanosoma indicum"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Indien goat fish"] <- "Parupeneus indicus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Indien ocean long nose pattor fish"] <- "Hipposcarus harid"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Indo Pacific Seargent"] <- "Abudefduf vaigiensis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Indo-Pacific"] <- "other"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Indo-pacific bonefish"] <- "Albula glossodonta"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Insular halfbeack"] <- "Hyporhamphus affinis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Insular halfbeak"] <- "Hyporhamphus affinis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Insular Halfbeak"] <- "Hyporhamphus affinis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Lembe"] <- "Johnius dorsalis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Lionfish spp."] <- "Pterois"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Lionfish Spp."] <- "Pterois"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Lionfish ssp."] <- "Pterois"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Logy sasa"] <- "other"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Long barbel goat fish"] <- "Parupeneus macronemus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Long Barbel Goat Fish"] <- "Parupeneus macronemus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Long Barbel Goatfish"] <- "Parupeneus macronemus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Long barrel goat fish"] <- "Parupeneus macronemus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Long barrel goatfish"] <- "Parupeneus macronemus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Long face emperor"] <- "Lethrinus microdon"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Long Face Emperor"] <- "Lethrinus microdon"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Long fish rudder fish"] <- "Kyphosus vaigiensis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Longbabbel goat fish"] <- "Parupeneus macronemus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "longfin bannerfish"] <- "Heniochus acuminatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Longfin bannerfish"] <- "Heniochus acuminatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Longfin rudder fish"] <- "Kyphosus vaigiensis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "low fin rudder fish"] <- "Kyphosus vaigiensis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Low fin rudder fish"] <- "Kyphosus vaigiensis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Low fin rudderfish"] <- "Kyphosus vaigiensis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Lowfin Ruddderfish"] <- "Kyphosus vaigiensis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Lowfin rudderfish"] <- "Kyphosus vaigiensis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Lowfin Rudderfish"] <- "Kyphosus vaigiensis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Lyretail grouper"] <- "Variola louti"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Manyspined angelfish"] <- "Centropyge multispinis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Marbled electric ray"] <- "Torpedo sinuspersici"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Mask bamboo fish"] <- "Heniochus monoceros"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Masked banner fish"] <- "Heniochus monoceros"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Miked species"] <- "other"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Mirex"] <- "other"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Mix (fish)"] <- "other"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Mix (Fish)"] <- "other"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Mix (Varilava)"] <- "other"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Mix species"] <- "other"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "mixed species"] <- "other"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Mixed species"] <- "other"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Mixed species of cardinal fish"] <- "Apogonidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Mixed Species"] <- "other"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Mixed spp."] <- "other"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Mooring idol"] <- "Zanclus cornutus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Moorish idol"] <- "Zanclus cornutus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Moorish Idol"] <- "Zanclus cornutus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Moray"] <- "Muraenidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Moray eel spp."] <- "Muraenidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Moray ssp."] <- "Muraenidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Morish Idol"] <- "Zanclus cornutus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Moustache conger"] <- "Conger cinereus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Moustache Conger"] <- "Conger cinereus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Moustache frigger"] <- "Balistoides viridescens"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Moustache Trigger Fish"] <- "Balistoides viridescens"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "moustashe trigger fish"] <- "Balistoides viridescens"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Mullet spp."] <- "Mugilidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "NAN"] <- "other"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Narrowbarred King mackrl"] <- "Scomberomorus commerson"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Needlefish spp."] <- "Belonidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Nigh fin rudder fish"] <- "Kyphosus cinerascens"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Nurse shark"] <- "Ginglymostomatidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "One spot snapper"] <- "Lutjanus monostigma"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Onespot snapper"] <- "Lutjanus monostigma"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Onespot Snapper"] <- "Lutjanus monostigma"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Orange line cardinal"] <- "Taeniamia fucata"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Orange Line Cardinal"] <- "Taeniamia fucata"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Orange line cardinal fish"] <- "Taeniamia fucata"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Orange Line Cardinal Fish"] <- "Taeniamia fucata"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Orange line cardinalfish"] <- "Taeniamia fucata"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Orange lined cardinal"] <- "Taeniamia fucata"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Orange line trigger"] <- "Balistapus undulatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Orange Lined Cardinal"] <- "Taeniamia fucata"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Orange lined cardinal fish"] <- "Taeniamia fucata"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Orange lined cardinalfish"] <- "Taeniamia fucata"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Orange spine unicornfish"] <- "Naso elegans"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Orange-lined cardinal"] <- "Taeniamia fucata"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Orange-lined cardinalfish"] <- "Taeniamia fucata"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Orange-lined Cardinalfish"] <- "Taeniamia fucata"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Orangelined Cardinal"] <- "Taeniamia fucata"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Palenose parroitfish"] <- "Scarus psittacus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Palenose parrot fish"] <- "Scarus psittacus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Picasso trigger fish"] <- "Rhinecanthus assasi"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Pink hindle barcuda"] <- "Sphyraena jello"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Pink-ear Emperor"] <- "Lethrinus lentjan"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Porcupinefish spp."] <- "Diodontidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Puffer fish"] <- "Tetraodontidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Puffer spp"] <- "Tetraodontidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Pufferfish spp."] <- "Tetraodontidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Ray (sol)"] <- "Batoidea"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Red chin flying fish"] <- "Exocoetidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "red gill emperor"] <- "Lethrinus rubrioperculatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Red gill emperor"] <- "Lethrinus rubrioperculatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Red Gill Emperor"] <- "Lethrinus rubrioperculatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Red rabbitfish"] <- "Siganidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Red spot goat fish"] <- "Parupeneus heptacanthus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Red spot goatfish"] <- "Parupeneus heptacanthus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Red Spot Goatfish"] <- "Parupeneus heptacanthus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Redchin flying fish"] <- "Exocoetidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Redchin flyingfish"] <- "Exocoetidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Redchin Flyingfish"] <- "Exocoetidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Redescent cardinalfish"] <- "Pristiapogon kallopterus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Redfin butterflyfish"] <- "Chaetodon trifasciatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Redfin Butterflyfish"] <- "Chaetodon trifasciatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Redgill Emperor"] <- "Lethrinus rubrioperculatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Redspot goatfish"] <- "Parupeneus heptacanthus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Redspot Goatfish"] <- "Parupeneus heptacanthus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Redspot goatfsih"] <- "Parupeneus heptacanthus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Reef  needle fish"] <- "Strongylura incisa"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Reef needlefish"] <- "Strongylura incisa"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Reef Needlefish"] <- "Strongylura incisa"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Reef neetle fish"] <- "Strongylura incisa"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Reef netdle fish"] <- "Strongylura incisa"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Reefneedle fish"] <- "Strongylura incisa"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Regal angel fish"] <- "Pygoplites diacanthus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Render spine mojarra"] <- "Gerreidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Ring tail cardinalfish"] <- "Ostorhinchus aureus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Ring tail surgeon fish"] <- "Acanthurus blochii"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Ring tail surgeonfish"] <- "Acanthurus blochii"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Ringed scorpionfish"] <- "Neomerinthe"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Ringtail surgeon fish"] <- "Acanthurus blochii"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Rosy goat fish"] <- "Parupeneus rubescens"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Rosy Goat fish"] <- "Parupeneus rubescens"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Rosy Goat Fish"] <- "Parupeneus rubescens"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Sabonto"] <- "Euthynnus affinis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Saddled parrotfish"] <- "Scarus scaber"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Sardine sp."] <- "Clupeidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "sardine spp."] <- "Clupeidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Sardine spp."] <- "Clupeidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Saringeba (mix)"] <- "other"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Sary geba"] <- "Glossogobius giuris"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Sary Geba"] <- "Glossogobius giuris"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Scissor tail pusilier"] <- "Caesio caerulaurea"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Scissor-tail Fusilier"] <- "Caesio caerulaurea"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Scissor-tail fusilierfish"] <- "Caesio caerulaurea"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Scissortail fusilier"] <- "Caesio caerulaurea"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Scissortail Fusilier"] <- "Caesio caerulaurea"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "scissortail fussilier"] <- "Caesio caerulaurea"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Scorpionfish spp."] <- "Neomerinthe"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Scorpionfish.spp"] <- "Neomerinthe"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Scribbled filefish"] <- "Aluterus scriptus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Scribbled snapper"] <- "Lutjanus rivulatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Scribbled Snapper"] <- "Lutjanus rivulatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Sea grasse parrot fish"] <- "Leptoscarus vaigiensis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Sea horse spp."] <- "Hippocampus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Seagrass parroitfish"] <- "Leptoscarus vaigiensis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "seagrass parrotfish"] <- "Leptoscarus vaigiensis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Seagrass parrotfish"] <- "Leptoscarus vaigiensis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Seagrass Parrotfish"] <- "Leptoscarus vaigiensis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Seagrasse parrot fish"] <- "Leptoscarus vaigiensis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Seagrasse Parrot Fish"] <- "Leptoscarus vaigiensis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Seagrasse parrotfish"] <- "Leptoscarus vaigiensis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Seagrasse Parrotfish"] <- "Leptoscarus vaigiensis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Semicercle angelfish"] <- "Pomacanthus semicirculatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Semicerlar angelfish"] <- "Pomacanthus semicirculatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Semicircle angel fish"] <- "Pomacanthus semicirculatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Semicircle Angel fish"] <- "Pomacanthus semicirculatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Shark spp."] <- "Elasmobranchii"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Side spot goat fish"] <- "Parupeneus pleurostigma"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Silver pompano"] <- "Trachinotus blochii"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Silver Pompano"] <- "Trachinotus blochii"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Slender prine mojarrah"] <- "Gerreidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Slender spine mojarra"] <- "Gerreidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Slender Spine Mojarra"] <- "Gerreidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Slender spine mosarra"] <- "Gerreidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Slenter spine majarra"] <- "Gerreidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Slinder spine mojarra"] <- "Gerreidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Slingfan wrasse"] <- "Epibulus insidiator"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Smooth tailled trevally"] <- "Selaroides leptolepis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Smooth-tailed trevally"] <- "Selaroides leptolepis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Smooth-tailed Trevally"] <- "Selaroides leptolepis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "smoothtail trevally"] <- "Selaroides leptolepis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Smoothtail trevally"] <- "Selaroides leptolepis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Smoth Mouth Fravailt"] <- "Carangidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Smoth tail trevally"] <- "Selaroides leptolepis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Snapper spp."] <- "Lutjanidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Snobe nose emperor"] <- "Lethrinus borbonicus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Snub nose emperor"] <- "Lethrinus borbonicus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Snubnosed Emperor"] <- "Lethrinus borbonicus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Soie"] <- "Soleidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Sole"] <- "Soleidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Sole spp."] <- "Soleidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Spadefish spp."] <- "Tripterodon orbis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Spine tooth parrot"] <- "Calotomus spinidens"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Spine tooth parrot fish"] <- "Calotomus spinidens"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Spinetooth parrotfish"] <- "Calotomus spinidens"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Spot fin lion fish"] <- "Pterois antennata"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Spotted butterflyfish"] <- "Chaetodon guttatissimus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "spotted Half Beak"] <- "Hemiramphus far"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Spotted half beak"] <- "Hemiramphus far"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Spotted sea hoose"] <- "Hippocampus kuda"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Spotted uni corn fish"] <- "Naso brevirostris"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Spotted unicorn fish"] <- "Naso brevirostris"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Spotted Unicorn Fish"] <- "Naso brevirostris"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Squirrel fish"] <- "Holocentridae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Squirrelfish spp."] <- "Holocentridae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Stallete rabbit fish"] <- "Siganus stellatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Star eye parroit fish"] <- "Calotomus carolinus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Star eye parrot fish"] <- "Calotomus carolinus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Star puffer"] <- "Arothron stellatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Stareye parroit"] <- "Calotomus carolinus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Stellate rubbitfish"] <- "Siganus stellatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Stellete rabbitfish"] <- "Siganus stellatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Stellete Rabbitfish"] <- "Siganus stellatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Stenler spine majorra"] <- "Gerreidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Sting Ray spp."] <- "Batoidea"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Stingray spp."] <- "Batoidea"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Stosky howkfish"] <- "Cirrhitus pinnulatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Stripe bristletooth"] <- "Ctenochaetus striatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Stripe clysinus"] <- "other"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Striped  bristletooth"] <- "Ctenochaetus striatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Striped brisile tooth"] <- "Ctenochaetus striatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Striped brisite"] <- "Ctenochaetus striatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Striped bristle tooth"] <- "Ctenochaetus striatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Striped Bristle Tooth"] <- "Ctenochaetus striatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Striped bristletooth"] <- "Ctenochaetus striatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Striped Bristletooth"] <- "Ctenochaetus striatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Striped cat fish"] <- "other"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Striped goatfish"] <- "Upeneus vittatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Tail spot squiddel fish"] <- "Sargocentron caudimaculatum"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Tail spot squirrel fish"] <- "Sargocentron caudimaculatum"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Tail spot squirrelfish"] <- "Sargocentron caudimaculatum"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "tailspot squirrelfish"] <- "Sargocentron caudimaculatum"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Tailspot squirrelfish"] <- "Sargocentron caudimaculatum"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Tailspot Squirrelfish"] <- "Sargocentron caudimaculatum"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Tailspot squirrelish"] <- "Sargocentron caudimaculatum"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Talang queen fish"] <- "Scomberoides commersonnianus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Talang Queen fish"] <- "Scomberoides commersonnianus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Thread fin buiter fly fish"] <- "Chaetodon auriga"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Thread fin butter fly fish"] <- "Chaetodon auriga"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Threadfin butter flyfish"] <- "Chaetodon auriga"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Thumbprint Spinecheek"] <- "Scolopsis affinis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Torovoky"] <- "Kyphosus cinerascens"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "trevally sp."] <- "Carangidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Trevally Sp."] <- "Carangidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "trevally spp"] <- "Carangidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Trevally spp"] <- "Carangidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "trevally spp."] <- "Carangidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Trevally spp."] <- "Carangidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Trevally Spp."] <- "Carangidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Trevelly spp"] <- "Carangidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Trevelly spp."] <- "Carangidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Tridescent cardinal fish"] <- "Pristiapogon kallopterus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Trigger fish"] <- "Balistidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Triggerfish"] <- "Balistidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Triggerfish spp."] <- "Balistidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "triggerfish ssp."] <- "Balistidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Triggerfish ssp."] <- "Balistidae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Triple Tail Whrasie"] <- "Cheilinus trilobatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Triple Tail Wrasse"] <- "Cheilinus trilobatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Tripletail wrabbe"] <- "Cheilinus trilobatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Tripltail wrass"] <- "Cheilinus trilobatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Tripple tail Wrasse"] <- "Cheilinus trilobatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Tripple Tail Wrasse"] <- "Cheilinus trilobatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Trippletail Wrasse"] <- "Cheilinus trilobatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Two barre goat fish"] <- "Parupeneus trifasciatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Two spot bristletooth"] <- "Ctenochaetus binotatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Two spot Bristletooth"] <- "Ctenochaetus binotatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Twobarred goatfish"] <- "Parupeneus trifasciatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Uanikoro sweeper"] <- "Pempheris vanicolensis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Unhulatet maray"] <- "Gymnothorax undulatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Unicornfish spp."] <- "Naso"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Unknow species"] <- "other"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Unknown spp."] <- "other"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Unsular Halfbeak"] <- "Hyporhamphus affinis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Vagabong Butter Fly Fish"] <- "Chaetodon vagabundus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Vagabong butter goat fish"] <- "Chaetodon vagabundus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Valala"] <- "other"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Vanikoro snapper"] <- "Pempheris vanicolensis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Varilava"] <- "Sicyopterus franouxi"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Varilava (mix)"] <- "Sicyopterus franouxi"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Varilava (Mix)"] <- "Sicyopterus franouxi"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Vary lava"] <- "Sicyopterus franouxi"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "vatsitsa"] <- "Chirocentrus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Vatsitsa"] <- "Chirocentrus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Wange el cat"] <- "other"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Wasse"] <- "Labridae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "While belly dansle"] <- "Amblyglyphidodon leucogaster"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "White fine wolf herrine"] <- "Chirocentrus nudus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "White margine unicorn fish"] <- "Naso annulatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "White margined anicornefish"] <- "Naso annulatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "White margined unicolor fish"] <- "Naso annulatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "White spotted guitar shark"] <- "Rhynchobatus djiddensis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "White Spotted Guitarfish"] <- "Rhynchobatus djiddensis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "White spotted gutarshark"] <- "Rhynchobatus djiddensis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "White spotted puffer"] <- "Arothron hispidus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "White Spotted Puffer"] <- "Arothron hispidus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "White-belly damsel"] <- "Amblyglyphidodon leucogaster"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "White-belly Damsel"] <- "Amblyglyphidodon leucogaster"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "White-spotted Guitarfish"] <- "Rhynchobatus djiddensis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Whitespotted puffer"] <- "Arothron hispidus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Whitespotted Puffer"] <- "Arothron hispidus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Whitespotted pufferfish"] <- "Arothron hispidus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Wire net file fish"] <- "Cantherhines pardalis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Wire net fire fish"] <- "Cantherhines pardalis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Wrasse"] <- "Labridae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Wrasse spp"] <- "Labridae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Wrasse spp."] <- "Labridae"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Yellow fin surgeon fish"] <- "Acanthurus xanthopterus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Yellow Fish Bream"] <- "Goldlined seabream"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Yellow Head Batter Fly Fish"] <- "Chaetodon xanthocephalus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Yellow lip emperor"] <- "Lethrinus xanthochilus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Yellow lips emperor"] <- "Lethrinus xanthochilus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Yellow margine trigger"] <- "Pseudobalistes flavimarginatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Yellow Margined Trigger"] <- "Pseudobalistes flavimarginatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Yellow spatted box"] <- "Ostracion cubicus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Yellow spotted emperor"] <- "Lethrinus erythracanthus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Yellow strife spotfish"] <- "Mulloidichthys flavolineatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Yellow Striped Goatfish"] <- "Mulloidichthys flavolineatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Yellow-Spotted trevally"] <- "Selaroides leptolepis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Yellow-spotted Trevally spp."] <- "Selaroides leptolepis"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Yellowmargined Mory"] <- "Gymnothorax flavimarginatus"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Yellowtop fusilier"] <- "Caesio xanthonota"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Zebra moray"] <- "Gymnomuraena zebra"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Zig Zag wrass"] <- "Halichoeres scapularis"

# Delete invertebrates and instances of no catch
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "Atoka tedro",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "Squid spp.",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "Bean turd sea cucumber",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "Carniloris Longis",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "Crab",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "cuttlefish",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "Cuttlefish",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "Cuttlefish spp.",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "Dorilisy",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "dorolisy",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "Dorolisy",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "Falaligalaky",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "Falalijaky",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "Golden Sandfish",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "Holothrua atra",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "holothuria atra",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[str_which(catch_mdg_tiliouine$species, "Holothuria", negate = TRUE),]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "Holothurie scabra",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "Lobster",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "Lobster spp.",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "Longlelegged spiny lobster",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "Maro mony (shell)",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "Murex",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[str_which(catch_mdg_tiliouine$species, "Octopus", negate = TRUE),]
catch_mdg_tiliouine <- catch_mdg_tiliouine[str_which(catch_mdg_tiliouine$species, "octopus", negate = TRUE),]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "Ornate long spine lobster",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "No (Ambalily)",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "No (fish)",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "No (Fish)",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "No (Octopus)",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "No (squid)",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "No (Squid)",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "no catch",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "No catch",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "No Catch",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "non",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "Non",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "none",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "None",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "NONE",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "Sea cucumber",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "Sea Cucumber",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "Sea cucumber spp.",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "Sea Cucumber spp.",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "Seacucumber",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "Seacucumber ssp.",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "Seacucumbers Acgia",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "Shell",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "Sicyopterus franouxi",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "Slipper lobster",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "Slipper lobster spp.",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "Slipper Lobster spp.",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "snail",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "Snail",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "squid",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "Squid",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "Squid Sp.",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "Squid Spp.",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "Steeper lobsto",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "Tsonkena",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "Zimele Papa",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "Zimely",]
catch_mdg_tiliouine <- catch_mdg_tiliouine[catch_mdg_tiliouine$species != "Zomely papa",]

# Find missepllings with taxize package
species_tiliouine <- unique(catch_mdg_tiliouine$species)
temp <- gnr_resolve(species_tiliouine, best_match_only = TRUE, canonical = TRUE)
temp <- filter(temp, temp$user_supplied_name != temp$matched_name2)

# Replace misspelled names
for(i in 1:nrow(temp)){
  catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == temp$user_supplied_name[i]] <- temp$matched_name2[i]
}

# Get WoRMS matches
species_tiliouine <- unique(catch_mdg_tiliouine$species)
worms_tiliouine <- worms_taxa2(species_tiliouine)

# Remove non-chordata
temp <- worms_tiliouine[worms_tiliouine$phylum != "Chordata", ]
catch_mdg_tiliouine <- filter(catch_mdg_tiliouine, !(catch_mdg_tiliouine$species %in% temp$supplied_name))
worms_tiliouine <- filter(worms_tiliouine, phylum == "Chordata")

# Remove sea turtles
worms_tiliouine <- filter(worms_tiliouine, family != "Cheloniidae")
catch_mdg_tiliouine <- filter(catch_mdg_tiliouine, species != "Chelonia mydas")

# Remove Poecilopsettidae because there is only one fish and the family has no nutrient estimates
worms_tiliouine <- filter(worms_tiliouine, valid_name != "Poecilopsettidae")
catch_mdg_tiliouine <- filter(catch_mdg_tiliouine, species != "Poecilopsettidae")

# Change "Pycnochromis leucurus" to "Pycnochromis" because the species is
# missing from the estimate table
worms_tiliouine$rank[worms_tiliouine$valid_name == "Pycnochromis leucurus"] <- "Genus"
worms_tiliouine$valid_name[worms_tiliouine$valid_name == "Pycnochromis leucurus"] <- "Pycnochromis"

# Change "Stegastes lacrymatus" to "Stegastes" because the species is missing
# from the estimate table
worms_tiliouine$rank[worms_tiliouine$valid_name == "Stegastes lacrymatus"] <- "Genus"
worms_tiliouine$valid_name[worms_tiliouine$valid_name == "Stegastes lacrymatus"] <- "Stegastes"
catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == "Stegastes lacrymatus"] <- "Stegastes"

# Correct names with WoRMS
temp <- filter(worms_tiliouine, worms_tiliouine$supplied_name != worms_tiliouine$valid_name)
for(i in 1:nrow(temp)){
  catch_mdg_tiliouine$species[catch_mdg_tiliouine$species == temp$supplied_name[i]] <- temp$valid_name[i]
}


## Weight column

# Delete no catch
catch_mdg_tiliouine <- filter(catch_mdg_tiliouine, weight_kg > 0)


## Length column

# Turn 0's to NA's
catch_mdg_tiliouine$length_cm[catch_mdg_tiliouine$length_cm == 0] <- NA

# Delete entries with length >= 200 (there are three and all are incorrect)
catch_mdg_tiliouine <- filter(catch_mdg_tiliouine, length_cm < 200)


## Gear column

# Make all entries lower case
catch_mdg_tiliouine$gear <- str_to_lower(catch_mdg_tiliouine$gear)

# Beach seine
catch_mdg_tiliouine$gear[catch_mdg_tiliouine$gear == "beach seine"] <- "beachseine"

# Boat seine
catch_mdg_tiliouine$gear[catch_mdg_tiliouine$gear == "boat sand"] <- "boatseine"
catch_mdg_tiliouine$gear[catch_mdg_tiliouine$gear == "boat sein"] <- "boatseine"
catch_mdg_tiliouine$gear[catch_mdg_tiliouine$gear == "boat seine"] <- "boatseine"

# Gill net
catch_mdg_tiliouine$gear[catch_mdg_tiliouine$gear == "gill net"] <- "gillnet"
catch_mdg_tiliouine$gear[catch_mdg_tiliouine$gear == "gill-net"] <- "gillnet"

# Spear
catch_mdg_tiliouine$gear[catch_mdg_tiliouine$gear == "harpon"] <- "spear"
catch_mdg_tiliouine$gear[catch_mdg_tiliouine$gear == "harpoon"] <- "spear"
catch_mdg_tiliouine$gear[catch_mdg_tiliouine$gear == "speagun"] <- "spear"
catch_mdg_tiliouine$gear[catch_mdg_tiliouine$gear == "spear gun"] <- "spear"
catch_mdg_tiliouine$gear[catch_mdg_tiliouine$gear == "speargon"] <- "spear"
catch_mdg_tiliouine$gear[catch_mdg_tiliouine$gear == "speargun"] <- "spear"
catch_mdg_tiliouine$gear[catch_mdg_tiliouine$gear == "spearn gun"] <- "spear"
catch_mdg_tiliouine$gear[catch_mdg_tiliouine$gear == "speat gun"] <- "spear"
catch_mdg_tiliouine$gear[catch_mdg_tiliouine$gear == "spoa gun"] <- "spear"

# Handline
catch_mdg_tiliouine$gear[catch_mdg_tiliouine$gear == "hook line"] <- "handline"
catch_mdg_tiliouine$gear[catch_mdg_tiliouine$gear == "hook-line"] <- "handline"
catch_mdg_tiliouine$gear[catch_mdg_tiliouine$gear == "hook'line"] <- "handline"
catch_mdg_tiliouine$gear[catch_mdg_tiliouine$gear == "hookline"] <- "handline"
catch_mdg_tiliouine$gear[catch_mdg_tiliouine$gear == "line"] <- "handline"

# Mosquito net
catch_mdg_tiliouine$gear[catch_mdg_tiliouine$gear == "mosquito gill net"] <- "mosquito net"

# NA
catch_mdg_tiliouine$gear[catch_mdg_tiliouine$gear == "nan"] <- NA
catch_mdg_tiliouine$gear[catch_mdg_tiliouine$gear == "tsaramaso"] <- NA
catch_mdg_tiliouine$gear[catch_mdg_tiliouine$gear == "vatsotsa"] <- NA
catch_mdg_tiliouine$gear[catch_mdg_tiliouine$gear == "votandra"] <- NA


## Add trip id's

# Empty column for trip id
catch_mdg_tiliouine$trip <- NA

# Unique combinations of date and effort
temp <- unique(catch_mdg_tiliouine[c("date", "effort", "gear")])

# Add empty trip column to temp
temp$trip <- NA

# Fill in trip column of temp
for(i in 1:nrow(temp)){
  temp$trip[i] <- trips[i]
}

# Delete used trip id's from trips
trips <- trips[(nrow(temp) + 1) : length(trips)]

# Fill in trip column of catch_mdg_tiliouine
for(i in 1:nrow(temp)){
  
  # If a gear is listed
  if(is.na(temp$gear[i]) == FALSE){
    
    # Row numbers that match this combination of date, effort, and gear
    x <- which(catch_mdg_tiliouine$date == temp$date[i] &
      catch_mdg_tiliouine$effort == temp$effort[i] &
      catch_mdg_tiliouine$gear == temp$gear[i])
    
  } else{ # If gear is NA
    
    # Row numbers that match this combination of date and effort
    x <- which(catch_mdg_tiliouine$date == temp$date[i] &
        catch_mdg_tiliouine$effort == temp$effort[i])
    
  }
  
  
  # Save trip id
  catch_mdg_tiliouine$trip[x] <- temp$trip[i]
  
}


## Add nutrient weights
catch_mdg_tiliouine <- nutrient_weights(catch_mdg_tiliouine, key = species_data, worms = worms_tiliouine)


## Add the L/Lopt indicator
catch_mdg_tiliouine <- lopt_indicator(catch_mdg_tiliouine, species_data, worms_tiliouine)


## Create trip data

# Add additional columns
catch_mdg_tiliouine$management <- NA
catch_mdg_tiliouine$country <- "madagascar"
catch_mdg_tiliouine$source <- "tiliouine2019"

# Create trip data
trip_tiliouine <- fish2trip(catch_mdg_tiliouine, dri = dri)

# Make date a date column
trip_tiliouine$date <- as.Date(trip_tiliouine$date, origin = "1970-01-01")


## Add to catch_mdg

# Aggregate by landing site and species
temp <- aggregate(weight_kg ~ species + site,
  data = catch_mdg_tiliouine, FUN = sum)

# Make a summary of landings by species by site
temp <- temp %>%
  pivot_wider(names_from = site, values_from = weight_kg)

# Bind to catch_mdg master
catch_mdg <- full_join(catch_mdg, temp)


## Add to trip data
trip <- bind_rows(trip, trip_tiliouine)
trip <- trip %>% select(-effort)


## Add to fish data
# fish <- bind_rows(fish, catch_mdg_tiliouine)
# fish <- fish %>% select(-(crew_id:zinc_mg))




##### Silva 2015 #####

## This section cleans and compiles data from:

## da Silva, I. M. (2015). Fisheries co-management: Ecological and social impacts. 
##    A case study of Northern Mozambique. Universidade de Aveiro.

## NB: There is no measurement of effort, so CPUE and nutrient yields will not be 
## available. I will add this to trip data with trip id's, but in this case the trip
## really just represents one sampling month (day?) for one gear at one location.
## this will allow reasonable estimates of nutrient concentrations.

## Import data
catch_moz_silva <- read_delim("data/raw-data/catch-monitoring/mozambique_catch-monitoring_silva2015.csv", 
    delim = ";", escape_double = FALSE, trim_ws = TRUE)

# Rename columns
catch_moz_silva <- catch_moz_silva %>%
  rename(site = estrato,
    length_cm = size)

# Select columns
catch_moz_silva <- catch_moz_silva %>%
  select(year, month, site, gear, species, length_cm, numsize)


## Date column

# Add an empty date column
catch_moz_silva$date <- NA

# Fill in date column
for(i in 1:nrow(catch_moz_silva)){
  
  catch_moz_silva$date[i] <- paste(catch_moz_silva$year[i], 
    catch_moz_silva$month[i],
    "01",
    sep = "-")
  
}

# Save as date
catch_moz_silva$date <- as.Date(catch_moz_silva$date)

# Remove year and month columns
catch_moz_silva <- catch_moz_silva %>%
  select(date, site, gear, species, length_cm, numsize)


## Site column

# Filter sites
catch_moz_silva <- catch_moz_silva %>% 
  filter(str_detect(catch_moz_silva$site, "ES", negate = TRUE)) %>% # remove estuarine sites
  filter(site != "Palma Sede Estuario") # remove Palma Sede estuarine site

# Correct site names
catch_moz_silva$site <- catch_moz_silva$site %>%
  str_remove_all(" MA") %>%
  str_remove_all(" 1") %>%
  str_remove_all(" 2") %>%
  str_remove_all(" Praia") %>%
  str_to_lower() %>%
  str_replace_all(pattern = "\xe9", replacement = "e")


## Gear column
catch_moz_silva$gear <- catch_moz_silva$gear %>%
  str_replace_all(pattern = "EMA", replacement = "other") %>%
  str_replace_all(pattern = "GAM", replacement = "other") %>%
  str_remove_all("1") %>%
  str_remove_all("2") %>%
  str_replace_all(pattern = "LME", replacement = "other") %>%
  str_replace_all(pattern = "MES", replacement = "other") %>%
  str_remove_all("gun") %>%
  str_replace_all(pattern = "line", replacement = "handline") %>%
  str_replace_all(pattern = "seine net", replacement = "beachseine") %>%
  str_replace_all(pattern = "Salema", replacement = "other") %>%
  str_to_lower()


## Species column

# Make sentence case and remove "sp." abbreviations
catch_moz_silva$species <- catch_moz_silva$species %>%
  str_to_sentence() %>%
  str_remove_all(" sp.")

# List of species
species_silva <- unique(catch_moz_silva$species)

# Find misspellings with taxize
temp <- gnr_resolve(species_silva, best_match_only = TRUE, canonical = TRUE)
temp <- filter(temp, temp$user_supplied_name != temp$matched_name2)

# Replace misspellings
for(i in 1:nrow(temp)){
  
  catch_moz_silva$species[catch_moz_silva$species == temp$user_supplied_name[i]] <-
    temp$matched_name2[i]
  
}

# Replace additional species by hand
catch_moz_silva$species <- catch_moz_silva$species %>%
  str_replace_all(pattern = "Calotomusnidens", replacement = "Calotomus spinidens") %>%
  str_replace_all(pattern = "Mossambica anguilla", replacement = "Anguilla mossambica")
  
# Remove unresolvable species
catch_moz_silva <- catch_moz_silva %>%
  filter(species != "Galeardo cuvier") %>%
  filter(species != "Outros") %>%
  filter(species != "Larvas")

# New list of species
species_silva <- unique(catch_moz_silva$species)

# Identify species in worms
worms_silva <- worms_taxa2(species_silva)

# Correct species names with worms
temp <- filter(worms_silva, worms_silva$supplied_name != worms_silva$valid_name)
for(i in 1:nrow(temp)){
  catch_moz_silva$species[catch_moz_silva$species == temp$supplied_name[i]] <-
    temp$valid_name[i]
}

# Remove non-chordata
temp <- filter(worms_silva, phylum != "Chordata")
for(i in 1:nrow(temp)){
  catch_moz_silva <- filter(catch_moz_silva, species != temp$valid_name[i])
}
worms_silva <- filter(worms_silva, phylum == "Chordata")


## Trip ID column

# Empty trip column
catch_moz_silva$trip <- NA

# Data frame of unique combinations of date, site, and gear
temp <- unique(catch_moz_silva[c("date", "site", "gear")])

# Fill in trip column
for(i in 1:nrow(temp)){
  
  catch_moz_silva$trip[catch_moz_silva$date == temp$date[i] &
      catch_moz_silva$site == temp$site[i] &
      catch_moz_silva$gear == temp$gear[i]] <- trips[i]
  
}

# Delete used trip id's
trips <- trips[(nrow(temp) + 1) : length(trips)]


## Mutate data to incorporate numsize column as new rows

# Replace NAs with 1
catch_moz_silva$numsize[which(is.na(catch_moz_silva$numsize))] <- 1

# Add additional rows
for(i in 1:nrow(catch_moz_silva)){
  
  # Number of observations of this species at this size
  x <- catch_moz_silva$numsize[i]
  
  # If there is more than one observation and new rows are needed
  if(x > 1){
    
    # Add additional rows
    catch_moz_silva <- rbind(catch_moz_silva, catch_moz_silva[rep(i, (x-1)), ])
    
  }
  
}

# Delete numsize column
catch_moz_silva <- select(catch_moz_silva, trip, date, site, gear, species, length_cm)



## Remove unrealistically big fishes

# Delete rows with na for length
catch_moz_silva <- catch_moz_silva[-which(is.na(catch_moz_silva$length_cm)), ]

# Empty vector to hold bad row numbers
temp <- c()

# Fill in temp
for(i in 1:nrow(worms_silva)){
  
  # If this is a species
  if(worms_silva$rank[i] == "Species"){
    
    # Extract max length
    ml <- species_data$max.length.tl[species_data$species == worms_silva$valid_name[i]]
    
  }
  
  # If this is a genus
  if(worms_silva$rank[i] == "Genus"){
    
    # Extract max length
    ml <- mean(species_data$max.length.tl[species_data$genus == worms_silva$valid_name[i]])
    
  }
  
  # If this is a family
  if(worms_silva$rank[i] == "Family"){
    
    # Extract max length
    ml <- mean(species_data$max.length.tl[species_data$family == worms_silva$valid_name[i]])
    
  }
  
  # Bad row numbers
  x <- which(catch_moz_silva$species == worms_silva$valid_name[i] &
      catch_moz_silva$length_cm > ml)
  
  # Save bad row numbers
  temp <- append(temp, x)
  
}

# Delete bad rows
catch_moz_silva <- catch_moz_silva[-temp, ]


## Add fish weights based on length and species

# Add empty weight column
catch_moz_silva$weight_kg <- NA

# Fill in weight column
catch_moz_silva <- fish_weight(catch_moz_silva, species_data, worms_silva)


## Add nutrient weights
catch_moz_silva <- nutrient_weights(catch_moz_silva, species_data, worms_silva)


## Add LLopt indicator
catch_moz_silva <- lopt_indicator(catch_moz_silva, species_data, worms_silva)


## Create trip data

# Add additional columns
catch_moz_silva$effort <- NA
catch_moz_silva$management <- NA
catch_moz_silva$country <- "mozambique"
catch_moz_silva$source<- "silva2015"

# Create trip data
trip_silva <- fish2trip(catch_moz_silva, dri = dri)

# Correct date
trip_silva$date <- as.Date(trip_silva$date, origin = "1970-01-01")


## Add to catch_moz

# Aggregate by landing site and species
temp <- aggregate(weight_kg ~ species + site,
  data = catch_moz_silva, FUN = sum)

# Make a summary of landings by species by site
catch_moz <- temp %>%
  pivot_wider(names_from = site, values_from = weight_kg)


## Add to trip data
trip <- bind_rows(trip, trip_silva)
trip <- trip %>% select(-effort)


## Add to fish data
fish <- bind_rows(fish, catch_moz_silva)
fish <- fish %>% select(-(calcium_mg:effort))




##### Exploration and Cleaning #####

## This section of the code explores the compiled trip data

## Date

# Sampling effort over time
ggplot(trip, aes(x = date, y = country)) +
  geom_point(alpha = 0.1) +
  theme_light()
ggplot(trip, aes(x = date, y = site)) +
  geom_point(alpha = 0.1) +
  theme_light()

# Sampling effort over gear
ggplot(trip, aes(x = date, y = gear)) +
  geom_point(alpha = 0.1) +
  theme_light()

# Landings over time
ggplot(trip, aes(x = date, y = cpue_kg.effort.day)) +
  geom_point(alpha = 0.1) +
  geom_smooth() +
  theme_light()

# Catch indicators over time
ggplot(trip, aes(x = date, y = LLopt)) +
  geom_point(alpha = 0.1) +
  geom_smooth() +
  theme_light()
ggplot(trip, aes(x = date, y = species_richness)) +
  geom_point(alpha = 0.1) +
  geom_smooth() +
  theme_light()

# Nutrients over time
ggplot(trip, aes(x = date, y = nutrients_pdv)) +
  geom_point(alpha = 0.1) +
  geom_smooth() +
  theme_light()
ggplot(trip, aes(x = date, y = nutrients_pdv.pue)) +
  geom_point(alpha = 0.1) +
  geom_smooth() +
  theme_light()


## Country

# Catch indicators by country
ggplot(trip, aes(x = country, y = cpue_kg.effort.day)) +
  geom_boxplot() +
  theme_light() +
  coord_cartesian(ylim = c(0, 10))
ggplot(trip, aes(x = country, y = LLopt)) +
  geom_boxplot() +
  theme_light() +
  coord_cartesian(ylim = c(0, 2))
ggplot(trip, aes(x = country, y = species_richness)) +
  geom_boxplot() +
  theme_light()

# Nutrients by country
ggplot(trip, aes(x = country, y = nutrients_pdv)) +
  geom_boxplot() +
  theme_light()
ggplot(trip, aes(x = country, y = nutrients_pdv.pue)) +
  geom_boxplot() +
  theme_light() +
  coord_cartesian(ylim = c(0, 100))
ggplot(trip, aes(x = country, y = calcium_mg.per.100g)) +
  geom_boxplot() +
  theme_light()
ggplot(trip, aes(x = country, y = calcium_mg.pue)) +
  geom_boxplot() +
  theme_light() +
  coord_cartesian(ylim = c(0, 10000))
ggplot(trip, aes(x = country, y = iron_mg.per.100g)) +
  geom_boxplot() +
  theme_light()
ggplot(trip, aes(x = country, y = iron_mg.pue)) +
  geom_boxplot() +
  theme_light() +
  coord_cartesian(ylim = c(0, 100))
ggplot(trip, aes(x = country, y = omega.3_g.per.100g)) +
  geom_boxplot() +
  theme_light()
ggplot(trip, aes(x = country, y = omega.3_g.pue)) +
  geom_boxplot() +
  theme_light() +
  coord_cartesian(ylim = c(0, 50))
ggplot(trip, aes(x = country, y = selenium_ug.per.100g)) +
  geom_boxplot() +
  theme_light()
ggplot(trip, aes(x = country, y = selenium_ug.pue)) +
  geom_boxplot() +
  theme_light() +
  coord_cartesian(ylim = c(0, 5000))
ggplot(trip, aes(x = country, y = vitamin.a_ug.per.100g)) +
  geom_boxplot() +
  theme_light()
ggplot(trip, aes(x = country, y = vitamin.a_ug.pue)) +
  geom_boxplot() +
  theme_light() +
  coord_cartesian(ylim = c(0, 10000))
ggplot(trip, aes(x = country, y = zinc_mg.per.100g)) +
  geom_boxplot() +
  theme_light()
ggplot(trip, aes(x = country, y = zinc_mg.pue)) +
  geom_boxplot() +
  theme_light() +
  coord_cartesian(ylim = c(0, 100))


## Site

# Catch indicators by site
ggplot(trip, aes(x = site, y = cpue_kg.effort.day)) +
  geom_boxplot(outlier.shape = NA) +
  theme_light() +
  coord_cartesian(ylim = c(0, 10))
ggplot(trip, aes(x = site, y = LLopt)) +
  geom_boxplot(outlier.shape = NA) +
  theme_light() +
  coord_cartesian(ylim = c(0, 2))

# Nutrients by site
ggplot(trip, aes(x = site, y = nutrients_pdv)) +
  geom_boxplot() +
  theme_light()
ggplot(trip, aes(x = site, y = nutrients_pdv.pue)) +
  geom_boxplot() +
  theme_light() +
  coord_cartesian(ylim = c(0, 100))
ggplot(trip, aes(x = site, y = calcium_mg.per.100g)) +
  geom_boxplot() +
  theme_light()
ggplot(trip, aes(x = site, y = calcium_mg.pue)) +
  geom_boxplot() +
  theme_light() +
  coord_cartesian(ylim = c(0, 10000))
ggplot(trip, aes(x = site, y = iron_mg.per.100g)) +
  geom_boxplot() +
  theme_light()
ggplot(trip, aes(x = site, y = iron_mg.pue)) +
  geom_boxplot() +
  theme_light() +
  coord_cartesian(ylim = c(0, 100))
ggplot(trip, aes(x = site, y = omega.3_g.per.100g)) +
  geom_boxplot() +
  theme_light()
ggplot(trip, aes(x = site, y = omega.3_g.pue)) +
  geom_boxplot() +
  theme_light() +
  coord_cartesian(ylim = c(0, 50))
ggplot(trip, aes(x = site, y = selenium_ug.per.100g)) +
  geom_boxplot() +
  theme_light()
ggplot(trip, aes(x = site, y = selenium_ug.pue)) +
  geom_boxplot() +
  theme_light() +
  coord_cartesian(ylim = c(0, 5000))
ggplot(trip, aes(x = site, y = vitamin.a_ug.per.100g)) +
  geom_boxplot() +
  theme_light()
ggplot(trip, aes(x = site, y = vitamin.a_ug.pue)) +
  geom_boxplot() +
  theme_light() +
  coord_cartesian(ylim = c(0, 10000))
ggplot(trip, aes(x = site, y = zinc_mg.per.100g)) +
  geom_boxplot() +
  theme_light()
ggplot(trip, aes(x = site, y = zinc_mg.pue)) +
  geom_boxplot() +
  theme_light() +
  coord_cartesian(ylim = c(0, 100))


## Gear

# Catch indicators by gear
ggplot(trip, aes(x = gear, y = cpue_kg.effort.day)) +
  geom_boxplot(outlier.shape = NA) +
  theme_light() +
  coord_cartesian(ylim = c(0, 10))
ggplot(trip, aes(x = gear, y = LLopt)) +
  geom_boxplot(outlier.shape = NA) +
  theme_light() +
  coord_cartesian(ylim = c(0, 2))
ggplot(fish, aes(x = gear, y = length_cm)) +
  geom_boxplot(outlier.shape = NA) +
  theme_light()

# Nutrients by gear
ggplot(trip, aes(x = gear, y = nutrients_pdv)) +
  geom_boxplot() +
  theme_light()
ggplot(trip, aes(x = gear, y = nutrients_pdv.pue)) +
  geom_boxplot() +
  theme_light() +
  coord_cartesian(ylim = c(0, 100))
ggplot(trip, aes(x = gear, y = calcium_mg.per.100g)) +
  geom_boxplot() +
  theme_light()
ggplot(trip, aes(x = gear, y = calcium_mg.pue)) +
  geom_boxplot() +
  theme_light() +
  coord_cartesian(ylim = c(0, 10000))
ggplot(trip, aes(x = gear, y = iron_mg.per.100g)) +
  geom_boxplot() +
  theme_light()
ggplot(trip, aes(x = gear, y = iron_mg.pue)) +
  geom_boxplot() +
  theme_light() +
  coord_cartesian(ylim = c(0, 100))
ggplot(trip, aes(x = gear, y = omega.3_g.per.100g)) +
  geom_boxplot() +
  theme_light()
ggplot(trip, aes(x = gear, y = omega.3_g.pue)) +
  geom_boxplot() +
  theme_light() +
  coord_cartesian(ylim = c(0, 50))
ggplot(trip, aes(x = gear, y = selenium_ug.per.100g)) +
  geom_boxplot() +
  theme_light()
ggplot(trip, aes(x = gear, y = selenium_ug.pue)) +
  geom_boxplot() +
  theme_light() +
  coord_cartesian(ylim = c(0, 5000))
ggplot(trip, aes(x = gear, y = vitamin.a_ug.per.100g)) +
  geom_boxplot() +
  theme_light()
ggplot(trip, aes(x = gear, y = vitamin.a_ug.pue)) +
  geom_boxplot() +
  theme_light() +
  coord_cartesian(ylim = c(0, 10000))
ggplot(trip, aes(x = gear, y = zinc_mg.per.100g)) +
  geom_boxplot() +
  theme_light()
ggplot(trip, aes(x = gear, y = zinc_mg.pue)) +
  geom_boxplot() +
  theme_light() +
  coord_cartesian(ylim = c(0, 100))


# Is CPUE related to nutrient concentrations?
ggplot(trip, aes(x = cpue_kg.effort.day, y = nutrients_pdv)) +
  geom_point(alpha = 0.3) +
  geom_smooth() +
  coord_cartesian(xlim = c(0, 10))
ggplot(trip, aes(x = cpue_kg.effort.day, y = calcium_mg.per.100g)) +
  geom_point(alpha = 0.3) +
  geom_smooth() +
  coord_cartesian(xlim = c(0, 10))
ggplot(trip, aes(x = cpue_kg.effort.day, y = iron_mg.per.100g)) +
  geom_point(alpha = 0.3) +
  geom_smooth() +
  coord_cartesian(xlim = c(0, 10))
ggplot(trip, aes(x = cpue_kg.effort.day, y = omega.3_g.per.100g)) +
  geom_point(alpha = 0.3) +
  geom_smooth() +
  coord_cartesian(xlim = c(0, 10))
ggplot(trip, aes(x = cpue_kg.effort.day, y = selenium_ug.per.100g)) +
  geom_point(alpha = 0.3) +
  geom_smooth() +
  coord_cartesian(xlim = c(0, 10))
ggplot(trip, aes(x = cpue_kg.effort.day, y = vitamin.a_ug.per.100g)) +
  geom_point(alpha = 0.3) +
  geom_smooth() +
  coord_cartesian(xlim = c(0, 10))
ggplot(trip, aes(x = cpue_kg.effort.day, y = zinc_mg.per.100g)) +
  geom_point(alpha = 0.3) +
  geom_smooth() +
  coord_cartesian(xlim = c(0, 10))


# Is LLopt related to nutrient concentrations?
ggplot(trip, aes(x = LLopt, y = nutrients_pdv)) +
  geom_point(alpha = 0.3) +
  geom_smooth() +
  coord_cartesian(xlim = c(0, 2))
ggplot(trip, aes(x = LLopt, y = calcium_mg.per.100g)) +
  geom_point(alpha = 0.3) +
  geom_smooth() +
  coord_cartesian(xlim = c(0, 2))
ggplot(trip, aes(x = LLopt, y = iron_mg.per.100g)) +
  geom_point(alpha = 0.3) +
  geom_smooth() +
  coord_cartesian(xlim = c(0, 2))
ggplot(trip, aes(x = LLopt, y = omega.3_g.per.100g)) +
  geom_point(alpha = 0.3) +
  geom_smooth() +
  coord_cartesian(xlim = c(0, 2))
ggplot(trip, aes(x = LLopt, y = selenium_ug.per.100g)) +
  geom_point(alpha = 0.3) +
  geom_smooth() +
  coord_cartesian(xlim = c(0, 2))
ggplot(trip, aes(x = LLopt, y = vitamin.a_ug.per.100g)) +
  geom_point(alpha = 0.3) +
  geom_smooth() +
  coord_cartesian(xlim = c(0, 2))
ggplot(trip, aes(x = LLopt, y = zinc_mg.per.100g)) +
  geom_point(alpha = 0.3) +
  geom_smooth() +
  coord_cartesian(xlim = c(0, 2))


# Is species richness related to nutrient concentrations?
ggplot(trip, aes(x = species_richness, y = nutrients_pdv)) +
  geom_point(alpha = 0.3) +
  geom_smooth()
ggplot(trip, aes(x = species_richness, y = calcium_mg.per.100g)) +
  geom_point(alpha = 0.3) +
  geom_smooth()
ggplot(trip, aes(x = species_richness, y = iron_mg.per.100g)) +
  geom_point(alpha = 0.3) +
  geom_smooth()
ggplot(trip, aes(x = species_richness, y = omega.3_g.per.100g)) +
  geom_point(alpha = 0.3) +
  geom_smooth()
ggplot(trip, aes(x = species_richness, y = selenium_ug.per.100g)) +
  geom_point(alpha = 0.3) +
  geom_smooth()
ggplot(trip, aes(x = species_richness, y = vitamin.a_ug.per.100g)) +
  geom_point(alpha = 0.3) +
  geom_smooth()
ggplot(trip, aes(x = species_richness, y = zinc_pdv)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm")


# LLopt ~ CPUE
ggplot(trip, aes(x = cpue_kg.effort.day, y = LLopt)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "lm") +
  theme_light() +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 4))
summary(lm(LLopt ~ cpue_kg.effort.day, trip))


# Species richness ~ LLopt
ggplot(trip, aes(x = LLopt, y = species_richness)) +
  geom_point(alpha = 0.1) +
  geom_smooth() +
  theme_light() +
  coord_cartesian(xlim = c(0, 4))

# Are some species overfished?

# 2/3 of observations and 50 most frequent species
fish2 <- filter(fish, species %in% names(sort(table(fish$species), decreasing = TRUE)[1:50]))

# Order species as a factor based on LLopt
fish2$species <- fct_reorder(fish2$species, fish2$LLopt)

# Some species are definitely overfished!
ggplot(fish2, aes(x = species, y = LLopt)) +
  geom_boxplot(outlier.alpha = 0.1) +
  geom_hline(yintercept = 1, color = "red") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# A closer look at the 15 most frequent species (still over half of observations)
fish3 <- filter(fish, species %in% names(sort(table(fish$species), decreasing = TRUE)[1:15]))
fish3$species <- fct_reorder(fish3$species, fish3$LLopt)
ggplot(fish3, aes(x = species, y = LLopt)) +
  geom_boxplot(outlier.alpha = 0.1) +
  geom_hline(yintercept = 1, color = "red") +
  ylab(expression(paste("Length : Optimum Length ", bgroup("(", frac(L, L[opt]), ")")))) +
  xlab("") +
  ggtitle("Growth overfishing in Kenya and Mozambique") +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "italic"),
    plot.margin = unit(c(1, 0.5, 0, 0.5), "cm"))
ggplot(fish3, aes(x = species, y = LLmat)) +
  geom_boxplot(outlier.alpha = 0.1) +
  geom_hline(yintercept = 1, color = "red") +
  ylab(expression(paste("Length : Length at Maturity ", bgroup("(", frac(L, L[mat]), ")")))) +
  xlab("") +
  ggtitle("Recruitment overfishing in Kenya and Mozambique") +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "italic"),
    plot.margin = unit(c(1, 0.5, 0, 0.5), "cm"))

# Nutrient and overall yields
ggplot(trip, aes(x = cpue_kg.effort.day, y = calcium_pdv.pue)) +
  geom_smooth(method = "lm", aes(color = "Calcium")) +
  geom_smooth(method = "lm",
    aes(x = cpue_kg.effort.day, y = iron_pdv.pue, color = "Iron")) +
  geom_smooth(method = "lm",
    aes(x = cpue_kg.effort.day, y = omega.3_pdv.pue, color = "Omega 3")) +
  geom_smooth(method = "lm",
    aes(x = cpue_kg.effort.day, y = selenium_pdv.pue, color = "Selenium")) +
  geom_smooth(method = "lm",
    aes(x = cpue_kg.effort.day, y = vitamin.a_pdv.pue, color = "Vitamin A")) +
  geom_smooth(method = "lm",
    aes(x = cpue_kg.effort.day, y = zinc_pdv.pue, color = "Zinc")) +
  scale_color_manual(name = "Nutrients",
    breaks=c('Calcium', 'Iron', 'Omega 3', 'Selenium', 'Vitamin A', 'Zinc'),
    values=c('Calcium'='orange', 'Iron'='skyblue', 'Omega 3'='pink', 'Selenium'='lightgreen',
      'Vitamin A'='purple', 'Zinc'='yellow')) +
  ylab(expression(paste("Nutrient Yield ", bgroup("(", frac("servings", "fisher x day"), ")")))) +
  xlab(expression(paste("Catch Per Unit Effort ", bgroup("(", frac(kg, "fisher x day"), ")")))) +
  theme_pubr() +
  theme(legend.background = element_blank()) +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 100))



##### Save Results #####

## Save Kenya catch summary

# Replace NAs with 0s
catch_kenya[is.na(catch_kenya)] <- 0

# Get mean values
catch_kenya <- catch_kenya %>%
  mutate(mean = rowMeans(across(kipini:waa), na.rm = TRUE))

# Make a proportional mean
catch_kenya$proportion <- catch_kenya$mean / sum(catch_kenya$mean, na.rm = TRUE)

# Save
write.csv(catch_kenya, "data/temp-data/03_CatchSummary_Kenya.csv", row.names = FALSE)


##  Save Madagascar catch summary

# Replace NAs with 0s
catch_mdg[is.na(catch_mdg)] <- 0

# Get mean values
catch_mdg <- catch_mdg %>%
  mutate(mean = rowMeans(across(andavadoaka:belavenoky), na.rm = TRUE))

# Make a proportional mean
catch_mdg$proportion <- catch_mdg$mean / sum(catch_kenya$mean, na.rm = TRUE)

# Save
write.csv(catch_mdg, "data/temp-data/03_CatchSummary_Madagascar.csv", row.names = FALSE)


## Save Mozambique catch summary

# Replace NAs with 0s
catch_moz[is.na(catch_moz)] <- 0

# Get mean values
catch_moz <- catch_moz %>%
  mutate(mean = rowMeans(across(guludo:vamizi), na.rm = TRUE))

# Make a proportional mean
catch_moz$proportion <- catch_moz$mean / sum(catch_moz$mean, na.rm = TRUE)

# Save
write.csv(catch_moz, "data/temp-data/03_CatchSummary_Mozambique.csv", row.names = FALSE)


## Save fish data
write.csv(fish, "data/clean-data/03_FishData.csv", row.names = FALSE)


## Save trip data

# Remove management column
trip <- select(trip, -"management")

# Make low frequency gears "other"
trip$gear[trip$gear == "net"] <- "other"
trip$gear[trip$gear == "nets"] <- "other"
trip$gear <- trip$gear %>%
  str_replace_all("castnet", "other") %>%
  str_replace_all("hookandstick", "other") %>%
  str_replace_all("mosquito net", "other") %>%
  str_replace_all("net-chachacha", "other") %>%
  str_replace_all("reefseine", "other") %>%
  str_replace_all("ringnet", "other") %>%
  str_replace_all("sardine net", "other") %>%
  str_replace_all("scoop net", "other") %>%
  str_replace_all("scoopnet", "other") %>%
  str_replace_all("trawling", "other") %>%
  str_replace_all("boatseine", "other") %>%
  str_replace_all("longline", "other")

# Save
write.csv(trip, "data/clean-data/03_TripData.csv", row.names = FALSE)

