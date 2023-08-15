## This script collates and cleans catch data to estimate the species composition of
##    artisanal catches by EEZ. It also explores the data and generates summary statistics.




##### Functions #####

## A function that consolidates Sea Around Us data to give a summary of landings (kg)
## by taxon. It consolidates taxa so that the least precise row subsumes all the more
## precise rows contained within it. For example, if one row is for a taxonomic class,
## any families, genera, or species within that class will be deleted and their
## landings added to the landings for that class. This is done on the presumption that
## members of said families, genera, or species could also be contained within the
## class observation.

## This function is only used within this workflow for countries for which we have
## catch monitoring data. Catch monitoring data will be used to add greater taxonomic
## precision to the consolidated SAU data, so there is a net gain in taxonomic precision.
## For countries without catch monitoring data, we will not consolidate SAU data.

## Parameter:
##    sau: a dataframe of 1-150 rows based on Sea Around Us data with landings aggregated by taxon
##          and unidentified and otherwise undesired taxa removed. Columns are
##          area_name, scientific_name, and tonnes.

## Output: a taxonomically consolidated version of the input.

consolidate_sau <- function(sau){
  
  # Get WORMS records for listed taxa
  if(nrow(sau) < 51){
    
    temp <- wm_records_taxamatch(sau$scientific_name)
    
    worms <- bind_rows(temp)
    
  }
  
  if(nrow(sau) > 50 && nrow(sau) < 101){
    
    temp <- wm_records_taxamatch(sau$scientific_name[1:50])
    temp2 <- wm_records_taxamatch(sau$scientific_name[51:nrow(sau)])
    
    worms <- bind_rows(c(temp, temp2))
    
  }
  
  if(nrow(sau) > 100){
    
    temp <- wm_records_taxamatch(sau$scientific_name[1:50])
    temp2 <- wm_records_taxamatch(sau$scientific_name[51:100])
    temp3 <- wm_records_taxamatch(sau$scientific_name[101:nrow(sau)])
    
    worms <- bind_rows(c(temp, temp2, temp3))
    
  }
  
  
  ## Delete secondary results for the taxa match (keeping only best matches)
  if (nrow(worms) > nrow(sau)){
    
    # Empty vector to hold list of row numbers to be deleted
    temp <- c()
    
    # Retrieve row numbers to be deleted
    for(i in 2:nrow(worms)){
      
      # test if this provided name matches the one in the row above it
      if (worms$scientificname[i] == worms$scientificname[i-1]){
        
        # if it is a duplicate, add to temp vector
        temp <- append(temp, i, after = length(temp))
      }
      
    }
    
    # Delete duplicate rows from worms
    worms <- worms[-temp, ]
    
  }
  
  # Combine WoRMS and SAU dataframes
  sau_revised <- bind_cols(worms, sau)
  
  # Send an error if the rows are not matched
  for (i in 1:nrow(sau_revised)){
    if (sau_revised$scientific_name[i] != sau_revised$scientificname[i]){
      stop(paste("Mismatched row in sau_revised (from SAU and WoRMS): ", i, sep = ""))
    }
  }
  
  # Remove all non chordata
  sau_revised <- filter(sau_revised, sau_revised$phylum == "Chordata")
  
  # Select only important columns of sau_revised
  sau_revised <- select(sau_revised, valid_AphiaID, valid_name, rank, kingdom, phylum, class,
    order, family, genus, tonnes)
  
  
  ## Consolidate partially resolved classes, orders, and families into the broadest classification
  
  # Temporary vectors to hold row numbers of taxa in need of consolidation
  classes <- c()
  orders <- c()
  families <- c()
  genera <- c()
  species <- c()
  
  # Identify which rows to consolidate to which taxonomic level
  for (i in 1:nrow(sau_revised)){
    
    # If the order is not known (must be consolidated to class)
    if (is.na(sau_revised$order[i]) == TRUE){
      classes <- append(classes, i, after = length(classes))
    }
    
    # If order is known but family is not (must be consolidated to order)
    if (is.na(sau_revised$order[i]) == FALSE && is.na(sau_revised$family[i]) == TRUE){
      orders <- append(orders, i, after = length(orders))
    }
    
    # If family is known but genus is not
    if (is.na(sau_revised$family[i]) == FALSE && is.na(sau_revised$genus[i]) == TRUE){
      families <- append(families, i, after = length(families))
    }
    
    # If genus is known but species is not
    if (sau_revised$rank[i] == "Genus"){
      genera <- append(genera, i, after = length(genera))
    }
    
    # If species is known
    if (sau_revised$rank[i] == "Species"){
      species <- append(species, i, after = length(species))
    }
    
  }
  
  
  ## Consolidate class
  
  # Vector of classes needing compilation
  classes <- sau_revised$class[classes]
  classes <- unique(classes)
  
  # Subset of sau_revised to be compiled by class
  sau_revised_class <- filter(sau_revised, class %in% classes)
  
  # Data frame with each class and total landings
  temp <- sau_revised_class %>%
    group_by(class) %>%
    summarize(landings = sum(tonnes))
  
  # Filter sau_revised_class for class level entries only
  sau_revised_class <- filter(sau_revised_class, rank == "Class")
  
  # Add correct landings to sau_revised_class
  for (i in 1:nrow(sau_revised_class)){
    
    # Get row number of this class in temp
    x <- match(sau_revised_class$class[i], temp$class)
    
    # Add landings to sau_revised_class
    sau_revised_class$tonnes[i] <- temp$landings[x]
  }
  
  
  ## Consolidate order
  
  # Vector of orders needing compilation
  orders <- sau_revised$order[orders]
  orders <- unique(orders)
  
  # Subset of sau_revised to be compiled by order
  sau_revised_order <- filter(sau_revised, order %in% orders) # Order level estimates fron SAU
  sau_revised_order <- filter(sau_revised_order, !(order %in% sau_revised_class$order)) # Remove if already included in class estimates
  
  # Data frame with each order and total landings
  temp <- sau_revised_order %>%
    group_by(order) %>%
    summarize(landings = sum(tonnes))
  
  # Filter sau_revised_order for order level entries only
  sau_revised_order <- filter(sau_revised_order, rank == "Order")
  
  # Add correct landings to sau_revised_order
  for (i in 1:nrow(sau_revised_order)){
    
    # Get row number of this class in temp
    x <- match(sau_revised_order$order[i], temp$order)
    
    # Add landings to sau_revised_order
    sau_revised_order$tonnes[i] <- temp$landings[x]
  }
  
  
  ## Consolidate families
  
  # Vector of families needing compilation
  families <- sau_revised$family[families]
  families <- unique(families)
  
  # Subset of sau_revised to be compiled by family
  sau_revised_family <- filter(sau_revised, family %in% families) # Family estimates
  sau_revised_family <- filter(sau_revised_family, !(order %in% sau_revised_order$order)) # not included in order estimates
  sau_revised_family <- filter(sau_revised_family, !(class %in% sau_revised_class$class)) # not included in class estimates
  
  # Data frame with each family and total landings
  temp <- sau_revised_family %>%
    group_by(family) %>%
    summarize(landings = sum(tonnes))
  
  # Filter family data frame for family level entries only
  sau_revised_family <- filter(sau_revised_family, rank == "Family")
  
  # Add correct landings to sau_revised_family
  for (i in 1:nrow(sau_revised_family)){
    
    # Get row number of this family in temp
    x <- match(sau_revised_family$family[i], temp$family)
    
    # Add landings to sau_revised_family
    sau_revised_family$tonnes[i] <- temp$landings[x]
  }
  
  
  ## Consolidate genera
  
  # Vector of genera needing compilation
  genera <- sau_revised$genus[genera]
  genera <- unique(genera)
  
  # Subset of sau_revised to be compiled by genus
  sau_revised_genus <- filter(sau_revised, genus %in% genera) # Genus estimates
  sau_revised_genus <- filter(sau_revised_genus, !(order %in% sau_revised_order$order)) # not included in order estimates
  sau_revised_genus <- filter(sau_revised_genus, !(class %in% sau_revised_class$class)) # not included in class estimates
  sau_revised_genus <- filter(sau_revised_genus, !(family %in% sau_revised_family$family)) # not included in family estimates
  
  # Data frame with each genus and total landings
  temp <- sau_revised_genus %>%
    group_by(genus) %>%
    summarize(landings = sum(tonnes))
  
  # Filter genus data frame for genus level entries only
  sau_revised_genus <- filter(sau_revised_genus, rank == "Genus")
  
  # Add correct landings to sau_revised_genus
  for (i in 1:nrow(sau_revised_genus)){
    
    # Get row number of this genus in temp
    x <- match(sau_revised_genus$genus[i], temp$genus)
    
    # Add landings to sau_revised_genus
    sau_revised_genus$tonnes[i] <- temp$landings[x]
  }
  
  
  ## Consolidate species
  
  # Vector of species needing compilation
  species <- sau_revised$valid_name[species]
  species <- unique(species)
  
  # Subset of sau_revised to be compiled by species
  sau_revised_species <- filter(sau_revised, valid_name %in% species) # Species estimates
  sau_revised_species <- filter(sau_revised_species, !(order %in% sau_revised_order$order)) # not included in order estimates
  sau_revised_species <- filter(sau_revised_species, !(class %in% sau_revised_class$class)) # not included in class estimates
  sau_revised_species <- filter(sau_revised_species, !(family %in% sau_revised_family$family)) # not included in family estimates
  sau_revised_species <- filter(sau_revised_species, !(genus %in% sau_revised_genus$genus)) # not included in genus estimates
  
  
  ## Combine consolidated estimates
  sau_revised <- bind_rows(sau_revised_class, sau_revised_order, sau_revised_family, sau_revised_genus, sau_revised_species)

  return(sau_revised)
  
}




## A function that combines Sea Around Us data with consolidated catch monitoring data to
## estimate a more precise taxonomic breakdown than what is available from SAU.

## Parameters are:
##  sau: a consolidated version of the Sea Around Us data, so that:
##        (1) there are no taxonomic overlaps (i.e., if there is a class listed in one row,
##        the landings from all families, genera, and species are added to the landings 
##        from that class and their rows are deleted). There are slight variations in the 
##        way SAU data are consolidated, so there is no function for that process. Instead, 
##        it is done de novo for each country in this script.
##        (2) the landings column is replaced with a ratio column, with each value being a
##        decimal estimate of the proportion of the overall catch accounted for by a taxon.
##        The column should sum to 1.
##        (3) the rank column is an ordered factor of taxonomic rank from Species < Class
##        (4) column names are: valid_AphiaID (from WoRMS database), valid_name, rank,
##        kingdom, phylum, class, order, family, genus, ratio
##  catch: a consolidated version of catch monitoring data with the same column names and
##        properties as sau, but with a landings column (kg) instead of ratio.

## Output is a dataframe that uses SAU catch composition estimates but improves it with
## taxonomic specificity from catch data, including relative catch composition of taxa
## within each SAU taxon.

improve_sau <- function(sau, catch){
  
  # Change column name for compatibility
  catch <- rename(catch, ratio = landings)
  
  ## Combine catch data with SAU

  # Empty data frame to contain new data
  temp <- sau[0,]
  
  for (i in 1:nrow(sau)){
    
    
    # If the SAU row is a taxonomic class
    if (sau$rank[i] == "Class"){
      
      # Subset catch data for members of this class
      temp2 <- filter(catch, class == sau$class[i])
      
      # If there are no matches, keep the row derived from SAU
      if (nrow(temp2) == 0){
        
        # Keep the SAU row
        temp <- bind_rows(temp, sau[i,])
        
      }
      
      # If there is only one match, choose the one that is more precise with SAU as default
      if (nrow(temp2) == 1){
        
        # If new row is more precise
        if (temp2$rank[1] < sau$rank[i]){
          
          # Assign correct ratio to new row
          temp2$ratio[1] <- sau$ratio[i]
          
          # Save the new row
          temp <- bind_rows(temp, temp2[1,])
          
        } else {
          
          # Keep the SAU row
          temp <- bind_rows(temp, sau[i,])
          
        }
        
      }
      
      # If there is more than one match
      if (nrow(temp2) > 1){
        
        # Turn landings into ratio
        temp2$ratio <- temp2$ratio / sum(temp2$ratio)
        
        # Calculate correct ratio as proportion of identified catch
        temp2$ratio <- temp2$ratio * sau$ratio[i]
        
        # Save the new row
        temp <- bind_rows(temp, temp2)
        
      }
      
    }
    
    
    # If the SAU row is a taxonomic order
    if (sau$rank[i] == "Order"){
      
      # Subset catch data for members of this order
      temp2 <- filter(catch, order == sau$order[i])
      
      # If there are no matches, keep the row derived from SAU
      if (nrow(temp2) == 0){
        
        # Keep the SAU row
        temp <- bind_rows(temp, sau[i,])
        
      }
      
      # If there is only one match, choose the one that is more precise with SAU as default
      if (nrow(temp2) == 1){
        
        # If new row is more precise
        if (temp2$rank[1] < sau$rank[i]){
          
          # Assign correct ratio to new row
          temp2$ratio[1] <- sau$ratio[i]
          
          # Save the new row
          temp <- bind_rows(temp, temp2[1,])
          
        } else {
          
          # Keep the SAU row
          temp <- bind_rows(temp, sau[i,])
          
        }
        
      }
      
      # If there is more than one match
      if (nrow(temp2) > 1){
        
        # Turn landings into ratio
        temp2$ratio <- temp2$ratio / sum(temp2$ratio)
        
        # Calculate correct ratio as proportion of identified catch
        temp2$ratio <- temp2$ratio * sau$ratio[i]
        
        # Save the new row
        temp <- bind_rows(temp, temp2)
        
      }
      
    }
    
    
    # If the SAU row is a taxonomic family
    if (sau$rank[i] == "Family"){
      
      # Subset catch data for members of this family
      temp2 <- filter(catch, family == sau$family[i])
      
      # If there are no matches, keep the row derived from SAU
      if (nrow(temp2) == 0){
        
        # Keep the SAU row
        temp <- bind_rows(temp, sau[i,])
        
      }
      
      # If there is only one match, choose the one that is more precise with SAU as default
      if (nrow(temp2) == 1){
        
        # If new row is more precise
        if (temp2$rank[1] < sau$rank[i]){
          
          # Assign correct ratio to new row
          temp2$ratio[1] <- sau$ratio[i]
          
          # Save the new row
          temp <- bind_rows(temp, temp2[1,])
          
        } else {
          
          # Keep the SAU row
          temp <- bind_rows(temp, sau[i,])
          
        }
        
      }
      
      # If there is more than one match
      if (nrow(temp2) > 1){
        
        # Turn landings into ratio
        temp2$ratio <- temp2$ratio / sum(temp2$ratio)
        
        # Calculate correct ratio as proportion of identified catch
        temp2$ratio <- temp2$ratio * sau$ratio[i]
        
        # Save the new row
        temp <- bind_rows(temp, temp2)
        
      }
      
    }
    
    
    # If the SAU row is a taxonomic genus
    if (sau$rank[i] == "Genus"){
      
      # Subset catch data for members of this genus
      temp2 <- filter(catch, genus == sau$genus[i])
      
      # If there are no matches, keep the row derived from SAU
      if (nrow(temp2) == 0){
        
        # Keep the SAU row
        temp <- bind_rows(temp, sau[i,])
        
      }
      
      # If there is only one match, choose the one that is more precise with SAU as default
      if (nrow(temp2) == 1){
        
        # If new row is more precise
        if (temp2$rank[1] < sau$rank[i]){
          
          # Assign correct ratio to new row
          temp2$ratio[1] <- sau$ratio[i]
          
          # Save the new row
          temp <- bind_rows(temp, temp2[1,])
          
        } else {
          
          # Keep the SAU row
          temp <- bind_rows(temp, sau[i,])
          
        }
        
      }
      
      # If there is more than one match
      if (nrow(temp2) > 1){
        
        # Turn landings into ratio
        temp2$ratio <- temp2$ratio / sum(temp2$ratio)
        
        # Calculate correct ratio as proportion of identified catch
        temp2$ratio <- temp2$ratio * sau$ratio[i]
        
        # Save the new rows
        temp <- bind_rows(temp, temp2)
        
      }
      
    }
    
    
    # If the SAU row is a species
    if (sau$rank[i] == "Species"){
      
      # Save the SAU row
      temp <- bind_rows(temp, sau[i,])
      
    }
    
  }
  
  # Remove AphiaID column
  temp <- select(temp, !valid_AphiaID)
  
  # Save national taxonomic breakdown for Tanzania
  return(temp)
    
}


## A function that pulls the WoRMS results you actually want: only the best matches
## for a list of taxa, returned as a dataframe with categories including Aphia ID,
## accepted name, supplied name, and full taxonomic information. This function can be
## used for lists of up to 500 species. NB: This is a wrapper for the
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
  
  if (length(species) > 350 && length(species) < 401){
    temp <- wm_records_taxamatch(species[1:50])
    temp2 <- wm_records_taxamatch(species[51:100])
    temp3 <- wm_records_taxamatch(species[101:150])
    temp4 <- wm_records_taxamatch(species[151:200])
    temp5 <- wm_records_taxamatch(species[201:250])
    temp6 <- wm_records_taxamatch(species[251:300])
    temp7 <- wm_records_taxamatch(species[301:350])
    temp8 <- wm_records_taxamatch(species[351:length(species)])
    df <- bind_rows(temp, temp2, temp3, temp4, temp5, temp6, temp7, temp8)
  }
  
  if (length(species) > 400 && length(species) < 451){
    temp <- wm_records_taxamatch(species[1:50])
    temp2 <- wm_records_taxamatch(species[51:100])
    temp3 <- wm_records_taxamatch(species[101:150])
    temp4 <- wm_records_taxamatch(species[151:200])
    temp5 <- wm_records_taxamatch(species[201:250])
    temp6 <- wm_records_taxamatch(species[251:300])
    temp7 <- wm_records_taxamatch(species[301:350])
    temp8 <- wm_records_taxamatch(species[351:400])
    temp9 <- wm_records_taxamatch(species[401:length(species)])
    df <- bind_rows(temp, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9)
  }
  
  if (length(species) > 450 && length(species) < 501){
    temp <- wm_records_taxamatch(species[1:50])
    temp2 <- wm_records_taxamatch(species[51:100])
    temp3 <- wm_records_taxamatch(species[101:150])
    temp4 <- wm_records_taxamatch(species[151:200])
    temp5 <- wm_records_taxamatch(species[201:250])
    temp6 <- wm_records_taxamatch(species[251:300])
    temp8 <- wm_records_taxamatch(species[351:400])
    temp9 <- wm_records_taxamatch(species[401:450])
    temp10 <- wm_records_taxamatch(species[451:length(species)])
    df <- bind_rows(temp, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10)
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




##### Kenya #####




##### * KE SAU #####

## The most recent Sea Around Us catch reconstruction for Kenya is:
##    McAlpine, A., & Zeller, D. (2020). Kenya: Updated catch reconstruction for 1950-2018.
##      In B. Derrick, M. Khalfallah, V. Relano, D. Zeller, & D. Pauly (Eds.), Updating to 2018
##      the 1950-2010 Marine Catch Reconstructions of the Sea Around Us: Part I – Africa,
##      Antarctica, Europe, and the North Atlantic (Vol. 5, pp. 46–59). Fisheries Center,
##      University of British Columbia. https://www.seaaroundus.org/data/#/eez/404

# Load SAU data for Kenya
sau_kenya <- read_csv("data/raw-data/sea-around-us/sau-kenya-11.25.22.csv")

# Filter for most recent year, SSF, and landings only (no discards)
sau_kenya <- filter(sau_kenya, year == max(sau_kenya$year))
sau_kenya <- filter(sau_kenya, catch_type == "Landings")
temp <- filter(sau_kenya, fishing_sector == "Artisanal")
temp2 <- filter(sau_kenya, fishing_sector == "Subsistence")
sau_kenya <- bind_rows(temp, temp2)

# Select relevant columns (country, scientific name, and landings)
sau_kenya <- select(sau_kenya, area_name, scientific_name, tonnes)

# Aggregate landings by taxon
sau_kenya <- aggregate(. ~ area_name + scientific_name, data = sau_kenya, sum)

# Remove unidentified, miscellaneous, and non-fish taxa
sau_kenya <- sau_kenya[-str_which(sau_kenya$scientific_name, "identified"), ]
sau_kenya <- sau_kenya[-str_which(sau_kenya$scientific_name, "Miscellaneous"), ]

# Consolidate sau_kenya
sau_kenya <- consolidate_sau(sau_kenya)

# Replace landings with ratios
sau_kenya$ratio <- sau_kenya$tonnes / sum(sau_kenya$tonnes)
sau_kenya <- select(sau_kenya, !tonnes)




##### * KE Catch Monitoring #####

# Load data
catch_kenya <- read_csv("data/temp-data/03_CatchSummary_Kenya.csv")

# Get worms results
worms_kenya <- worms_taxa2(catch_kenya$species)

# Combine worms with catch data
catch_kenya <- bind_cols(worms_kenya, catch_kenya)
catch_kenya <- select(catch_kenya, valid_name, rank:genus, mean)
catch_kenya <- rename(catch_kenya, landings = mean)

# Make rank columns an ordered factor so they can be compared
sau_kenya$rank <- factor(sau_kenya$rank, levels = c("Species", "Genus", "Family", "Order", "Class"), ordered = TRUE)
catch_kenya$rank <- factor(catch_kenya$rank, levels = c("Species", "Genus", "Family", "Order", "Class"), ordered = TRUE)

# Improve taxonomic precision of SAU data
kenya <- improve_sau(sau_kenya, catch_kenya)




##### Tanzania #####




##### * TZ SAU #####

## The most recent Sea Around Us catch reconstruction for Tanzania is:
##    White, R., Page, E., & Noël, S.-L. (2020). Tanzania: Updated catch reconstruction
##      for 2011-2018. In B. Derrick, M. Khalfallah, V. Relano, D. Zeller, & D. Pauly
##      (Eds.), Updating to 2018 the 1950-2010 marine catch reconstructions of the Sea 
##      Around Us: Part I – Africa, Antarctica, Europe, and the North Atlantic (Vol. 5, 
##      pp. 77–80). Fisheries Center, University of British Columbia. 
##      https://www.seaaroundus.org/data/#/eez/834


# Load SAU data for Tanzania
sau_tanzania <- read_csv("data/raw-data/sea-around-us/sau-tanzania-12.7.22.csv")

# Filter for most recent year, SSF, and landings only (no discards)
sau_tanzania <- filter(sau_tanzania, year == max(sau_tanzania$year))
sau_tanzania <- filter(sau_tanzania, catch_type == "Landings")
temp <- filter(sau_tanzania, fishing_sector == "Artisanal")
temp2 <- filter(sau_tanzania, fishing_sector == "Subsistence")
sau_tanzania <- bind_rows(temp, temp2)

# Select relevant columns (country, scientific name, and landings)
sau_tanzania <- select(sau_tanzania, area_name, scientific_name, tonnes)

# Aggregate landings by taxon
sau_tanzania <- aggregate(. ~ area_name + scientific_name, data = sau_tanzania, sum)

# Remove unidentified taxa
sau_tanzania <- sau_tanzania[-str_which(sau_tanzania$scientific_name, "identified"), ]

# Consolidate sau_tanzania
sau_tanzania <- consolidate_sau(sau_tanzania)

# Replace landings with ratios
sau_tanzania$ratio <- sau_tanzania$tonnes / sum(sau_tanzania$tonnes)
sau_tanzania <- select(sau_tanzania, !tonnes)




##### * TZ Catch Monitoring #####

## Import catch data for Tanzania
catch_tanzania <- read_csv("data/raw-data/catch-monitoring/tanzania-catch-monitoring-tzmlf2020.csv")
catch_tanzania <- select(catch_tanzania, Species, Total)

# Remove freshwater and unidentified taxa
catch_tanzania <- filter(catch_tanzania, Species != "Clarias gariepin")
catch_tanzania <- filter(catch_tanzania, Species != "Other")

# Get worms results
worms_tanzania <- worms_taxa2(catch_tanzania$Species)

# Combine worms with catch data
catch_tanzania <- bind_cols(worms_tanzania, catch_tanzania)
catch_tanzania <- select(catch_tanzania, valid_name:genus, Total)
catch_tanzania <- rename(catch_tanzania, landings = Total)

# Make rank columns an ordered factor so they can be compared
sau_tanzania$rank <- factor(sau_tanzania$rank, levels = c("Species", "Genus", "Family", "Order", "Class"), ordered = TRUE)
catch_tanzania$rank <- factor(catch_tanzania$rank, levels = c("Species", "Genus", "Family", "Order", "Class"), ordered = TRUE)

# Improve taxonomic precision of SAU data
tanzania <- improve_sau(sau_tanzania, catch_tanzania)




##### Mozambique #####




##### * MZ SAU #####

## The most recent Sea Around Us catch reconstruction for Mozambique is:
##    Vianna, G. M. S. (2020). Mozambique: Updated catch reconstruction for 2011–2018.
##        In B. Derrick, M. Khalfallah, V. Relano, D. Zeller, & D. Pauly (Eds.), Updating 
##        to 2018 the 1950-2010 marine catch reconstructions of the Sea Around Us: Part 
##        I – Africa, Antarctica, Europe, and the North Atlantic (Vol. 5, pp. 65–68). 
##        Fisheries Center, University of British Columbia. 
##        https://www.seaaroundus.org/data/#/eez/508



# Load SAU data for mozambique
sau_mozambique <- read_csv("data/raw-data/sea-around-us/sau-mozambique-12.7.22.csv")

# Filter for most recent year, SSF, and landings only (no discards)
sau_mozambique <- filter(sau_mozambique, year == max(sau_mozambique$year))
sau_mozambique <- filter(sau_mozambique, catch_type == "Landings")
temp <- filter(sau_mozambique, fishing_sector == "Artisanal")
temp2 <- filter(sau_mozambique, fishing_sector == "Subsistence")
sau_mozambique <- bind_rows(temp, temp2)

# Select relevant columns (country, scientific name, and landings)
sau_mozambique <- select(sau_mozambique, area_name, scientific_name, tonnes)

# Aggregate landings by taxon
sau_mozambique <- aggregate(. ~ area_name + scientific_name, data = sau_mozambique, sum)

# Remove unidentified taxa
sau_mozambique <- sau_mozambique[-str_which(sau_mozambique$scientific_name, "identified"), ]

# Consolidate sau_mozambique
sau_mozambique <- consolidate_sau(sau_mozambique)

# Replace landings with ratios
sau_mozambique$ratio <- sau_mozambique$tonnes / sum(sau_mozambique$tonnes)
sau_mozambique <- select(sau_mozambique, !tonnes)




##### * MZ Catch Monitoring #####

# Load data
catch_moz <- read_csv("data/temp-data/03_CatchSummary_Mozambique.csv")

# Get worms results
worms_moz <- worms_taxa2(catch_moz$species)

# Combine worms with catch data
catch_moz <- bind_cols(worms_moz, catch_moz)
catch_moz <- select(catch_moz, valid_name, rank:genus, mean)
catch_moz <- rename(catch_moz, landings = mean)

# Make rank columns an ordered factor so they can be compared
sau_mozambique$rank <- factor(sau_mozambique$rank, levels = c("Species", "Genus", "Family", "Order", "Class"), ordered = TRUE)
catch_moz$rank <- factor(catch_moz$rank, levels = c("Species", "Genus", "Family", "Order", "Class"), ordered = TRUE)

# Improve taxonomic precision of SAU data
moz <- improve_sau(sau_mozambique, catch_moz)

# Save national taxonomic breakdown for moz
mozambique <- moz




##### Madagascar #####




##### * MDG SAU #####

## The most recent Sea Around Us catch reconstruction for Madagascar is:
##    White, R., Noël, S.-L., Christ, H., Relano, V., Sicnawa, F., & Tsui, G. (2020). 
##        Madagascar and smaller islands in the Western Indian Ocean: Updated catch 
##        reconstructions for 2011-2018. In B. Derrick, M. Khalfallah, V. Relano, D. 
##        Zeller, & D. Pauly (Eds.), Updating to 2018 the 1950-2010 marine catch 
##        reconstructions of the Sea Around Us: Part I - Africa, Antarctica, Europe, 
##        and the North Atlantic. (Vol. 5, pp. 81–99). Fisheries Center, University of 
##        British Columbia. https://sau-technical-reports.s3.amazonaws.com/174_175_251_252_450_480_638_690_972_White_et_al_2020_Madagascar_and_surrounding_islands_FCRR.pdf?Signature=AgaeGD2Mv7hctEu%2BYI1IHlhAsQc%3D&Expires=1670592067&AWSAccessKeyId=ASIAWJJOTCLMPJHMKW7W&x-amz-security-token=IQoJb3JpZ2luX2VjEB0aCXVzLXdlc3QtMiJHMEUCIQCNK8cYTCbSyc%2B9AC9sxK44sxNovl9Kh/sgobjv3Jac%2BwIgbh7QlsTWrLPqDQ2eeCRNGo60X%2BhMt4J9zGFxwu6x2H0qzAQIRhADGgw0MzIyNzg4MDMxNjAiDHR5gIcrG23s3%2BJPcSqpBCTYfzrZCOagNotLXzfRhDiwgp7y77nJ8t9TqB1JFQZoK/kyBZI7SgseG1vBA5GBdT8oK5MM6/5nH%2BqrsWGoZqgHTqBdDd8JWGHpbYYjEj600RQMEFoUHouDCompkEccwS4l0jxM97PLBqX%2BB8BsLfykZRuep5Jwp2lKdk8lzONDQTSkpsNgM%2Bj5yHjPPfkOE0sxc68S4Ovcr4maEDyWheLBEdZ9HsmuzuIj5NHV8PX3MVfwJ7A4hyEIhQtCX9uUoJ3yNqMNV/14%2BWEaopxSY4bli04EUmCFIRjp2ctT9qGBGvhDziWDcagn35dS%2BSb/cckPaUSsLI/OOETKcgNszGWqCR%2BznbwrO5uB6DJ4bZNx7u7ZKrszCYIR6DSsUQdxZt0cppIMuIAb5TVk9ALRDMAdwOhHJE2DA%2BkSo9rHnhn2QmoaXIgiyzMuiPs/97Xx5OYQoxbxPz1y2B1mJ8DvF8xFL2nZ9YfHotsCj6dIsR93AaS76Sj3uzrmJofwPMQjNP6Yuv28l53e8Yrilp9%2BkD1D4ACvW9LQMhKZBa1Rft6BflQhesxEzZTzPgHZDB2PfuoMxR7HwRA%2B50jgL8EcAzEKN4WB7Es3k5lTu%2BWxUEsJo2WwTMY6Vd5HEN87ECmeWF%2Bn4Whl4/kPsKz0QSJRIr948V/c78NMotTOqfCz3%2B8ZVZZqF%2BP9X9cb5tF6uGyaJSQAUXroydJXMuIz%2BMJqDKfrCo%2Bp2fKe6F0wtMTHnAY6qQH6u5FifccNsnvtnehEgZOYn2G5NZr8DvVt/fSFpRZfLV5dWybGMMRJosLSy/ZSsV7hEwH76dVZmz9AiVbtngrYj4RMSFOPT23s1%2B8Z9h%2BQDzNhagxY29Kl8/KTshihz6JRJqE/0usT972uaGy30s7HYINnl6qsQBx3lQJ%2BJpvlCcdJFFLdczQBKx3dQseXM6GyDRQtDZBJMNijUWd1N5H2hvPho4gaZwdz

# Load SAU data for madagascar
sau_madagascar <- read_csv("data/raw-data/sea-around-us/sau-madagascar-12.8.22.csv")

# Filter for most recent year, SSF, and landings only (no discards)
sau_madagascar <- filter(sau_madagascar, year == max(sau_madagascar$year))
sau_madagascar <- filter(sau_madagascar, catch_type == "Landings")
temp <- filter(sau_madagascar, fishing_sector == "Artisanal")
temp2 <- filter(sau_madagascar, fishing_sector == "Subsistence")
sau_madagascar <- bind_rows(temp, temp2)

# Select relevant columns (country, scientific name, and landings)
sau_madagascar <- select(sau_madagascar, area_name, scientific_name, tonnes)

# Aggregate landings by taxon
sau_madagascar <- aggregate(. ~ area_name + scientific_name, data = sau_madagascar, sum)

# Remove unidentified taxa
sau_madagascar <- sau_madagascar[-str_which(sau_madagascar$scientific_name, "identified"), ]
sau_madagascar <- sau_madagascar[-str_which(sau_madagascar$scientific_name, "Miscellaneous"), ]

# Consolidate sau_madagascar
sau_madagascar <- consolidate_sau(sau_madagascar)

# Replace landings with ratios
sau_madagascar$ratio <- sau_madagascar$tonnes / sum(sau_madagascar$tonnes)
sau_madagascar <- select(sau_madagascar, !tonnes)




##### * MDG Catch Data #####

# Load data
catch_mdg <- read_csv("data/temp-data/03_CatchSummary_Madagascar.csv")

# Get worms results
worms_mdg <- worms_taxa2(catch_mdg$species)

# Delete unidentified species
catch_mdg <- catch_mdg[catch_mdg$species %in% worms_mdg$valid_name, ]

# Combine worms with catch data
catch_mdg <- bind_cols(worms_mdg, catch_mdg)
catch_mdg <- select(catch_mdg, valid_name, rank:genus, mean)
catch_mdg <- rename(catch_mdg, landings = mean)

# Make rank columns an ordered factor so they can be compared
sau_madagascar$rank <- factor(sau_madagascar$rank, levels = c("Species", "Genus", "Family", "Order", "Class"), ordered = TRUE)
catch_mdg$rank <- factor(catch_mdg$rank, levels = c("Species", "Genus", "Family", "Order", "Class"), ordered = TRUE)

# Improve taxonomic precision of SAU data
mdg <- improve_sau(sau_madagascar, catch_mdg)

# Save national taxonomic breakdown for mdg
madagascar <- mdg

# Remove Batoidae
madagascar$ratio[madagascar$valid_name == "Elasmobranchii"] <- madagascar$ratio[madagascar$valid_name == "Elasmobranchii"] + madagascar$ratio[madagascar$valid_name == "Batoidea"]
madagascar <- filter(madagascar, valid_name != "Batoidea")

# Make Ferdauia orthogrammus family level
madagascar$rank[madagascar$valid_name == "Ferdauia orthogrammus"] <- "Family"



##### Maldives #####




##### * MD SAU #####

## The most recent Sea Around Us catch reconstruction for Maldives is:
##    White, R., Hood, L., Derrick, B., Relano, V., & Zeller, D. (2020). South Asia and 
##        Indian Ocean island catch updates to 2018. In B. Derrick, M. Khalfallah, V. 
##        Relano, D. Zeller, & D. Pauly (Eds.), Updating to 2018 the 1950-2010 marine 
##        catch reconstructions of the Sea Around Us. Part II: The Americas and Asia-
##        Pacific (Vol. 6, pp. 365–382). Fisheries Center, University of British Columbia.
##        https://sau-technical-reports.s3.amazonaws.com/50_86_144_356_357_462_586_White_et_al_2020_South_Asian_Indian_Ocean_Islands_FCRR.pdf?Signature=SKj1dlKD%2BtG6aK4nI%2Bvtk%2BmLgrw%3D&Expires=1670592495&AWSAccessKeyId=ASIAWJJOTCLMPJHMKW7W&x-amz-security-token=IQoJb3JpZ2luX2VjEB0aCXVzLXdlc3QtMiJHMEUCIQCNK8cYTCbSyc%2B9AC9sxK44sxNovl9Kh/sgobjv3Jac%2BwIgbh7QlsTWrLPqDQ2eeCRNGo60X%2BhMt4J9zGFxwu6x2H0qzAQIRhADGgw0MzIyNzg4MDMxNjAiDHR5gIcrG23s3%2BJPcSqpBCTYfzrZCOagNotLXzfRhDiwgp7y77nJ8t9TqB1JFQZoK/kyBZI7SgseG1vBA5GBdT8oK5MM6/5nH%2BqrsWGoZqgHTqBdDd8JWGHpbYYjEj600RQMEFoUHouDCompkEccwS4l0jxM97PLBqX%2BB8BsLfykZRuep5Jwp2lKdk8lzONDQTSkpsNgM%2Bj5yHjPPfkOE0sxc68S4Ovcr4maEDyWheLBEdZ9HsmuzuIj5NHV8PX3MVfwJ7A4hyEIhQtCX9uUoJ3yNqMNV/14%2BWEaopxSY4bli04EUmCFIRjp2ctT9qGBGvhDziWDcagn35dS%2BSb/cckPaUSsLI/OOETKcgNszGWqCR%2BznbwrO5uB6DJ4bZNx7u7ZKrszCYIR6DSsUQdxZt0cppIMuIAb5TVk9ALRDMAdwOhHJE2DA%2BkSo9rHnhn2QmoaXIgiyzMuiPs/97Xx5OYQoxbxPz1y2B1mJ8DvF8xFL2nZ9YfHotsCj6dIsR93AaS76Sj3uzrmJofwPMQjNP6Yuv28l53e8Yrilp9%2BkD1D4ACvW9LQMhKZBa1Rft6BflQhesxEzZTzPgHZDB2PfuoMxR7HwRA%2B50jgL8EcAzEKN4WB7Es3k5lTu%2BWxUEsJo2WwTMY6Vd5HEN87ECmeWF%2Bn4Whl4/kPsKz0QSJRIr948V/c78NMotTOqfCz3%2B8ZVZZqF%2BP9X9cb5tF6uGyaJSQAUXroydJXMuIz%2BMJqDKfrCo%2Bp2fKe6F0wtMTHnAY6qQH6u5FifccNsnvtnehEgZOYn2G5NZr8DvVt/fSFpRZfLV5dWybGMMRJosLSy/ZSsV7hEwH76dVZmz9AiVbtngrYj4RMSFOPT23s1%2B8Z9h%2BQDzNhagxY29Kl8/KTshihz6JRJqE/0usT972uaGy30s7HYINnl6qsQBx3lQJ%2BJpvlCcdJFFLdczQBKx3dQseXM6GyDRQtDZBJMNijUWd1N5H2hvPho4gaZwdz


# Load SAU data for maldives
sau_maldives <- read_csv("data/raw-data/sea-around-us/sau-maldives-12.8.22.csv")

# Filter for most recent year, SSF, and landings only (no discards)
sau_maldives <- filter(sau_maldives, year == max(sau_maldives$year))
sau_maldives <- filter(sau_maldives, catch_type == "Landings")
temp <- filter(sau_maldives, fishing_sector == "Artisanal")
temp2 <- filter(sau_maldives, fishing_sector == "Subsistence")
sau_maldives <- bind_rows(temp, temp2)

# Select relevant columns (country, scientific name, and landings)
sau_maldives <- select(sau_maldives, area_name, scientific_name, tonnes)

# Aggregate landings by taxon
sau_maldives <- aggregate(. ~ area_name + scientific_name, data = sau_maldives, sum)

# Remove unidentified taxa
sau_maldives <- sau_maldives[-str_which(sau_maldives$scientific_name, "identified"), ]

# Replace landings with ratios
sau_maldives$ratio <- sau_maldives$tonnes / sum(sau_maldives$tonnes)
sau_maldives <- select(sau_maldives, !tonnes)

# Get worms results
worms_maldives <- worms_taxa2(sau_maldives$scientific_name)

# Combine worms with catch data
maldives <- bind_cols(worms_maldives, sau_maldives)
maldives <- select(maldives, valid_name:genus, ratio)

# Make Ferdauia orthogrammus family level
maldives$rank[maldives$valid_name == "Ferdauia orthogrammus"] <- "Family"


##### Seychelles #####




##### * SCH SAU #####

## The most recent catch reconstruction for Seychelles is:
##    Christ, H. J., White, R., Hood, L., Vianna, G. M. S., & Zeller, D. (2020). A 
##        baseline for the blue economy: Catch and effort history in the Republic of 
##        Seychelles’ domestic fisheries. Frontiers in Marine Science, 7. 
##        https://doi.org/10.3389/fmars.2020.00269
##    White, R., Noël, S.-L., Christ, H., Relano, V., Sicnawa, F., & Tsui, G. (2020). 
##        Madagascar and smaller islands in the Western Indian Ocean: Updated catch 
##        reconstructions for 2011-2018. In B. Derrick, M. Khalfallah, V. Relano, D. 
##        Zeller, & D. Pauly (Eds.), Updating to 2018 the 1950-2010 marine catch 
##        reconstructions of the Sea Around Us: Part I - Africa, Antarctica, Europe, 
##        and the North Atlantic. (Vol. 5, pp. 81–99). Fisheries Center, University of 
##        British Columbia. 
##        https://sau-technical-reports.s3.amazonaws.com/174_175_251_252_450_480_638_690_972_White_et_al_2020_Madagascar_and_surrounding_islands_FCRR.pdf?Signature=AgaeGD2Mv7hctEu%2BYI1IHlhAsQc%3D&Expires=1670592067&AWSAccessKeyId=ASIAWJJOTCLMPJHMKW7W&x-amz-security-token=IQoJb3JpZ2luX2VjEB0aCXVzLXdlc3QtMiJHMEUCIQCNK8cYTCbSyc%2B9AC9sxK44sxNovl9Kh/sgobjv3Jac%2BwIgbh7QlsTWrLPqDQ2eeCRNGo60X%2BhMt4J9zGFxwu6x2H0qzAQIRhADGgw0MzIyNzg4MDMxNjAiDHR5gIcrG23s3%2BJPcSqpBCTYfzrZCOagNotLXzfRhDiwgp7y77nJ8t9TqB1JFQZoK/kyBZI7SgseG1vBA5GBdT8oK5MM6/5nH%2BqrsWGoZqgHTqBdDd8JWGHpbYYjEj600RQMEFoUHouDCompkEccwS4l0jxM97PLBqX%2BB8BsLfykZRuep5Jwp2lKdk8lzONDQTSkpsNgM%2Bj5yHjPPfkOE0sxc68S4Ovcr4maEDyWheLBEdZ9HsmuzuIj5NHV8PX3MVfwJ7A4hyEIhQtCX9uUoJ3yNqMNV/14%2BWEaopxSY4bli04EUmCFIRjp2ctT9qGBGvhDziWDcagn35dS%2BSb/cckPaUSsLI/OOETKcgNszGWqCR%2BznbwrO5uB6DJ4bZNx7u7ZKrszCYIR6DSsUQdxZt0cppIMuIAb5TVk9ALRDMAdwOhHJE2DA%2BkSo9rHnhn2QmoaXIgiyzMuiPs/97Xx5OYQoxbxPz1y2B1mJ8DvF8xFL2nZ9YfHotsCj6dIsR93AaS76Sj3uzrmJofwPMQjNP6Yuv28l53e8Yrilp9%2BkD1D4ACvW9LQMhKZBa1Rft6BflQhesxEzZTzPgHZDB2PfuoMxR7HwRA%2B50jgL8EcAzEKN4WB7Es3k5lTu%2BWxUEsJo2WwTMY6Vd5HEN87ECmeWF%2Bn4Whl4/kPsKz0QSJRIr948V/c78NMotTOqfCz3%2B8ZVZZqF%2BP9X9cb5tF6uGyaJSQAUXroydJXMuIz%2BMJqDKfrCo%2Bp2fKe6F0wtMTHnAY6qQH6u5FifccNsnvtnehEgZOYn2G5NZr8DvVt/fSFpRZfLV5dWybGMMRJosLSy/ZSsV7hEwH76dVZmz9AiVbtngrYj4RMSFOPT23s1%2B8Z9h%2BQDzNhagxY29Kl8/KTshihz6JRJqE/0usT972uaGy30s7HYINnl6qsQBx3lQJ%2BJpvlCcdJFFLdczQBKx3dQseXM6GyDRQtDZBJMNijUWd1N5H2hvPho4gaZwdz


# Load SAU data for Seychelles
sau_seychelles <- read_csv("data/raw-data/sea-around-us/sau-seychelles-12.8.22.csv")

# Filter for most recent year, SSF, and landings only (no discards)
sau_seychelles <- filter(sau_seychelles, year == max(sau_seychelles$year))
sau_seychelles <- filter(sau_seychelles, catch_type == "Landings")
temp <- filter(sau_seychelles, fishing_sector == "Artisanal")
temp2 <- filter(sau_seychelles, fishing_sector == "Subsistence")
sau_seychelles <- bind_rows(temp, temp2)

# Select relevant columns (country, scientific name, and landings)
sau_seychelles <- select(sau_seychelles, area_name, scientific_name, tonnes)

# Aggregate landings by taxon
sau_seychelles <- aggregate(. ~ area_name + scientific_name, data = sau_seychelles, sum)

# Remove unidentified taxa
sau_seychelles <- sau_seychelles[-str_which(sau_seychelles$scientific_name, "identified"), ]

# Consolidate sau_seychelles
sau_seychelles <- consolidate_sau(sau_seychelles)

# Replace landings with ratios
sau_seychelles$ratio <- sau_seychelles$tonnes / sum(sau_seychelles$tonnes)
sau_seychelles <- select(sau_seychelles, !tonnes)




##### * SCH Catch Data #####

## Import catch data
catch_seychelles <- read_excel("data/raw-data/catch-monitoring/seychelles-catch-monitoring-sfa.xlsx", 
  col_types = c("numeric", "text", "text", 
    "numeric", "text", "skip"))
catch_seychelles <- select(catch_seychelles, Scientific_name, Total)

# Rename columns
catch_seychelles <- rename(catch_seychelles, c(scientific_name = Scientific_name, landings = Total))

# Correct Selachimorpha to Elasmobranchii
catch_seychelles$scientific_name <- replace(catch_seychelles$scientific_name,
  list = which(catch_seychelles$scientific_name == "Selachimorpha (Pleurotremata)"),
  values = "Elasmobranchii")

# Correct Thunnini to Scombridae
catch_seychelles$scientific_name <- replace(catch_seychelles$scientific_name,
  list = which(catch_seychelles$scientific_name == "Thunnini"),
  values = "Scombridae")

# Aggregate by taxon
catch_seychelles <- aggregate(. ~ scientific_name, data = catch_seychelles, sum)

# Get WORMS records for listed taxa
temp <- wm_records_taxamatch(catch_seychelles$scientific_name[1:50])
temp2 <- wm_records_taxamatch(catch_seychelles$scientific_name[51:100])
temp3 <- wm_records_taxamatch(catch_seychelles$scientific_name[101:nrow(catch_seychelles)])
worms_seychelles <- bind_rows(temp, temp2, temp3)

## Delete secondary results for the taxa match (keeping only best matches)
if (nrow(worms_seychelles) > nrow(catch_seychelles)){
  
  # Empty vector to hold list of row numbers to be deleted
  temp <- c()
  
  # Retrieve row numbers to be deleted
  for(i in 2:nrow(worms_seychelles)){
  
    # test if this provided name matches the one in the row above it
    if (worms_seychelles$scientificname[i] == worms_seychelles$scientificname[i-1]){
  
      # if it is a duplicate, add to temp vector
      temp <- append(temp, i, after = length(temp))
    }
  
  }
  
  # Delete duplicate rows from worms_seychelles
  worms_seychelles <- worms_seychelles[-temp, ]
  
}

# Select important columns of worms_seychelles
worms_seychelles <- select(worms_seychelles, valid_AphiaID, valid_name, rank,
  kingdom, phylum, class, order, family, genus)

# Combine worms with catch data
catch_seychelles <- bind_cols(worms_seychelles, catch_seychelles)
catch_seychelles <- select(catch_seychelles, !scientific_name)

# Make rank columns an ordered factor so they can be compared
sau_seychelles$rank <- factor(sau_seychelles$rank, levels = c("Species", "Genus", "Family", "Order", "Class"), ordered = TRUE)
catch_seychelles$rank <- factor(catch_seychelles$rank, levels = c("Species", "Genus", "Family", "Order", "Class"), ordered = TRUE)

# Improve taxonomic precision of SAU data
seychelles <- improve_sau(sau_seychelles, catch_seychelles)

# Make Percoidei an order
seychelles$rank[seychelles$valid_name == "Percoidei"] <- "Order"



##### Reunion #####




##### * RN SAU #####

## The most recent catch reconstruction for Reunion is:
##    White, R., Noël, S.-L., Christ, H., Relano, V., Sicnawa, F., & Tsui, G. (2020). 
##        Madagascar and smaller islands in the Western Indian Ocean: Updated catch 
##        reconstructions for 2011-2018. In B. Derrick, M. Khalfallah, V. Relano, D. 
##        Zeller, & D. Pauly (Eds.), Updating to 2018 the 1950-2010 marine catch 
##        reconstructions of the Sea Around Us: Part I - Africa, Antarctica, Europe, 
##        and the North Atlantic. (Vol. 5, pp. 81–99). Fisheries Center, University of 
##        British Columbia. 
##        https://sau-technical-reports.s3.amazonaws.com/174_175_251_252_450_480_638_690_972_White_et_al_2020_Madagascar_and_surrounding_islands_FCRR.pdf?Signature=AgaeGD2Mv7hctEu%2BYI1IHlhAsQc%3D&Expires=1670592067&AWSAccessKeyId=ASIAWJJOTCLMPJHMKW7W&x-amz-security-token=IQoJb3JpZ2luX2VjEB0aCXVzLXdlc3QtMiJHMEUCIQCNK8cYTCbSyc%2B9AC9sxK44sxNovl9Kh/sgobjv3Jac%2BwIgbh7QlsTWrLPqDQ2eeCRNGo60X%2BhMt4J9zGFxwu6x2H0qzAQIRhADGgw0MzIyNzg4MDMxNjAiDHR5gIcrG23s3%2BJPcSqpBCTYfzrZCOagNotLXzfRhDiwgp7y77nJ8t9TqB1JFQZoK/kyBZI7SgseG1vBA5GBdT8oK5MM6/5nH%2BqrsWGoZqgHTqBdDd8JWGHpbYYjEj600RQMEFoUHouDCompkEccwS4l0jxM97PLBqX%2BB8BsLfykZRuep5Jwp2lKdk8lzONDQTSkpsNgM%2Bj5yHjPPfkOE0sxc68S4Ovcr4maEDyWheLBEdZ9HsmuzuIj5NHV8PX3MVfwJ7A4hyEIhQtCX9uUoJ3yNqMNV/14%2BWEaopxSY4bli04EUmCFIRjp2ctT9qGBGvhDziWDcagn35dS%2BSb/cckPaUSsLI/OOETKcgNszGWqCR%2BznbwrO5uB6DJ4bZNx7u7ZKrszCYIR6DSsUQdxZt0cppIMuIAb5TVk9ALRDMAdwOhHJE2DA%2BkSo9rHnhn2QmoaXIgiyzMuiPs/97Xx5OYQoxbxPz1y2B1mJ8DvF8xFL2nZ9YfHotsCj6dIsR93AaS76Sj3uzrmJofwPMQjNP6Yuv28l53e8Yrilp9%2BkD1D4ACvW9LQMhKZBa1Rft6BflQhesxEzZTzPgHZDB2PfuoMxR7HwRA%2B50jgL8EcAzEKN4WB7Es3k5lTu%2BWxUEsJo2WwTMY6Vd5HEN87ECmeWF%2Bn4Whl4/kPsKz0QSJRIr948V/c78NMotTOqfCz3%2B8ZVZZqF%2BP9X9cb5tF6uGyaJSQAUXroydJXMuIz%2BMJqDKfrCo%2Bp2fKe6F0wtMTHnAY6qQH6u5FifccNsnvtnehEgZOYn2G5NZr8DvVt/fSFpRZfLV5dWybGMMRJosLSy/ZSsV7hEwH76dVZmz9AiVbtngrYj4RMSFOPT23s1%2B8Z9h%2BQDzNhagxY29Kl8/KTshihz6JRJqE/0usT972uaGy30s7HYINnl6qsQBx3lQJ%2BJpvlCcdJFFLdczQBKx3dQseXM6GyDRQtDZBJMNijUWd1N5H2hvPho4gaZwdz

# Load SAU data for Reunion
sau_reunion <- read_csv("data/raw-data/sea-around-us/sau-reunion-12.8.22.csv")

# Filter for most recent year, SSF, and landings only (no discards)
sau_reunion <- filter(sau_reunion, year == max(sau_reunion$year))
sau_reunion <- filter(sau_reunion, catch_type == "Landings")
temp <- filter(sau_reunion, fishing_sector == "Artisanal")
temp2 <- filter(sau_reunion, fishing_sector == "Subsistence")
sau_reunion <- bind_rows(temp, temp2)

# Select relevant columns (country, scientific name, and landings)
sau_reunion <- select(sau_reunion, area_name, scientific_name, tonnes)

# Aggregate landings by taxon
sau_reunion <- aggregate(. ~ area_name + scientific_name, data = sau_reunion, sum)

# Remove unidentified taxa
sau_reunion <- sau_reunion[-str_which(sau_reunion$scientific_name, "identified"), ]

# Consolidate sau_reunion
sau_reunion <- consolidate_sau(sau_reunion)

# Replace landings with ratios
sau_reunion$ratio <- sau_reunion$tonnes / sum(sau_reunion$tonnes)
sau_reunion <- select(sau_reunion, !tonnes)




##### * RN Catch Data #####

## Import catch data for Reunion
catch_reunion <- read_excel("data/raw-data/catch-monitoring/reunion-catch-monitoring-biais1992.xlsx")
catch_reunion <- select(catch_reunion, scientific_name, catch_kg)

# Get WORMS records for listed taxa
temp <- wm_records_taxamatch(catch_reunion$scientific_name[1:50])
temp2 <- wm_records_taxamatch(catch_reunion$scientific_name[51:100])
temp3 <- wm_records_taxamatch(catch_reunion$scientific_name[101:nrow(catch_reunion)])
worms_reunion <- bind_rows(temp, temp2, temp3)

## Delete secondary results for the taxa match (keeping only best matches)
if (nrow(worms_reunion) > nrow(catch_reunion)){
  
  # Empty vector to hold list of row numbers to be deleted
  temp <- c()
  
  # Retrieve row numbers to be deleted
  for(i in 2:nrow(worms_reunion)){
  
    # test if this provided name matches the one in the row above it
    if (worms_reunion$scientificname[i] == worms_reunion$scientificname[i-1]){
  
      # if it is a duplicate, add to temp vector
      temp <- append(temp, i, after = length(temp))
    }
  
  }
  
  # Delete duplicate rows from worms_reunion
  worms_reunion <- worms_reunion[-temp, ]
  
}

# Select important columns of worms_reunion
worms_reunion <- select(worms_reunion, valid_AphiaID, valid_name, rank,
  kingdom, phylum, class, order, family, genus)

# Combine worms with catch data
catch_reunion <- bind_cols(worms_reunion, catch_reunion)
catch_reunion <- select(catch_reunion, !scientific_name)
catch_reunion <- rename(catch_reunion, landings = catch_kg)

# Make rank columns an ordered factor so they can be compared
sau_reunion$rank <- factor(sau_reunion$rank, levels = c("Species", "Genus", "Family", "Order", "Class"), ordered = TRUE)
catch_reunion$rank <- factor(catch_reunion$rank, levels = c("Species", "Genus", "Family", "Order", "Class"), ordered = TRUE)

# Improve taxonomic precision of SAU data
reunion <- improve_sau(sau_reunion, catch_reunion)




##### Mayotte #####




##### * MY SAU #####

## The most recent catch reconstruction for Mayotte is:
##    White, R., Noël, S.-L., Christ, H., Relano, V., Sicnawa, F., & Tsui, G. (2020). 
##        Madagascar and smaller islands in the Western Indian Ocean: Updated catch 
##        reconstructions for 2011-2018. In B. Derrick, M. Khalfallah, V. Relano, D. 
##        Zeller, & D. Pauly (Eds.), Updating to 2018 the 1950-2010 marine catch 
##        reconstructions of the Sea Around Us: Part I - Africa, Antarctica, Europe, 
##        and the North Atlantic. (Vol. 5, pp. 81–99). Fisheries Center, University of 
##        British Columbia. https://sau-technical-reports.s3.amazonaws.com/174_175_251_252_450_480_638_690_972_White_et_al_2020_Madagascar_and_surrounding_islands_FCRR.pdf?Signature=AgaeGD2Mv7hctEu%2BYI1IHlhAsQc%3D&Expires=1670592067&AWSAccessKeyId=ASIAWJJOTCLMPJHMKW7W&x-amz-security-token=IQoJb3JpZ2luX2VjEB0aCXVzLXdlc3QtMiJHMEUCIQCNK8cYTCbSyc%2B9AC9sxK44sxNovl9Kh/sgobjv3Jac%2BwIgbh7QlsTWrLPqDQ2eeCRNGo60X%2BhMt4J9zGFxwu6x2H0qzAQIRhADGgw0MzIyNzg4MDMxNjAiDHR5gIcrG23s3%2BJPcSqpBCTYfzrZCOagNotLXzfRhDiwgp7y77nJ8t9TqB1JFQZoK/kyBZI7SgseG1vBA5GBdT8oK5MM6/5nH%2BqrsWGoZqgHTqBdDd8JWGHpbYYjEj600RQMEFoUHouDCompkEccwS4l0jxM97PLBqX%2BB8BsLfykZRuep5Jwp2lKdk8lzONDQTSkpsNgM%2Bj5yHjPPfkOE0sxc68S4Ovcr4maEDyWheLBEdZ9HsmuzuIj5NHV8PX3MVfwJ7A4hyEIhQtCX9uUoJ3yNqMNV/14%2BWEaopxSY4bli04EUmCFIRjp2ctT9qGBGvhDziWDcagn35dS%2BSb/cckPaUSsLI/OOETKcgNszGWqCR%2BznbwrO5uB6DJ4bZNx7u7ZKrszCYIR6DSsUQdxZt0cppIMuIAb5TVk9ALRDMAdwOhHJE2DA%2BkSo9rHnhn2QmoaXIgiyzMuiPs/97Xx5OYQoxbxPz1y2B1mJ8DvF8xFL2nZ9YfHotsCj6dIsR93AaS76Sj3uzrmJofwPMQjNP6Yuv28l53e8Yrilp9%2BkD1D4ACvW9LQMhKZBa1Rft6BflQhesxEzZTzPgHZDB2PfuoMxR7HwRA%2B50jgL8EcAzEKN4WB7Es3k5lTu%2BWxUEsJo2WwTMY6Vd5HEN87ECmeWF%2Bn4Whl4/kPsKz0QSJRIr948V/c78NMotTOqfCz3%2B8ZVZZqF%2BP9X9cb5tF6uGyaJSQAUXroydJXMuIz%2BMJqDKfrCo%2Bp2fKe6F0wtMTHnAY6qQH6u5FifccNsnvtnehEgZOYn2G5NZr8DvVt/fSFpRZfLV5dWybGMMRJosLSy/ZSsV7hEwH76dVZmz9AiVbtngrYj4RMSFOPT23s1%2B8Z9h%2BQDzNhagxY29Kl8/KTshihz6JRJqE/0usT972uaGy30s7HYINnl6qsQBx3lQJ%2BJpvlCcdJFFLdczQBKx3dQseXM6GyDRQtDZBJMNijUWd1N5H2hvPho4gaZwdz

# Load SAU data for Mayotte
sau_mayotte <- read_csv("data/raw-data/sea-around-us/sau-mayotte-12.8.22.csv")

# Filter for most recent year, SSF, and landings only (no discards)
sau_mayotte <- filter(sau_mayotte, year == max(sau_mayotte$year))
sau_mayotte <- filter(sau_mayotte, catch_type == "Landings")
temp <- filter(sau_mayotte, fishing_sector == "Artisanal")
temp2 <- filter(sau_mayotte, fishing_sector == "Subsistence")
sau_mayotte <- bind_rows(temp, temp2)

# Select relevant columns (country, scientific name, and landings)
sau_mayotte <- select(sau_mayotte, area_name, scientific_name, tonnes)

# Aggregate landings by taxon
sau_mayotte <- aggregate(. ~ area_name + scientific_name, data = sau_mayotte, sum)

# Remove unidentified taxa
sau_mayotte <- sau_mayotte[-str_which(sau_mayotte$scientific_name, "identified"), ]

# Remove Scylla serrata
sau_mayotte <- sau_mayotte[sau_mayotte$scientific_name != "Scylla serrata", ]

# Replace landings with ratios
sau_mayotte$ratio <- sau_mayotte$tonnes / sum(sau_mayotte$tonnes)
sau_mayotte <- select(sau_mayotte, !tonnes)

# Get worms results
worms_mayotte <- worms_taxa2(sau_mayotte$scientific_name)

# Combine worms with catch data
mayotte <- bind_cols(worms_mayotte, sau_mayotte)
mayotte <- select(mayotte, valid_name:genus, ratio)




##### Mauritius #####




##### * MR SAU #####

## The most recent catch reconstruction for Mauritius is:
##    White, R., Noël, S.-L., Christ, H., Relano, V., Sicnawa, F., & Tsui, G. (2020). 
##        Madagascar and smaller islands in the Western Indian Ocean: Updated catch 
##        reconstructions for 2011-2018. In B. Derrick, M. Khalfallah, V. Relano, D. 
##        Zeller, & D. Pauly (Eds.), Updating to 2018 the 1950-2010 marine catch 
##        reconstructions of the Sea Around Us: Part I - Africa, Antarctica, Europe, 
##        and the North Atlantic. (Vol. 5, pp. 81–99). Fisheries Center, University of 
##        British Columbia. https://sau-technical-reports.s3.amazonaws.com/174_175_251_252_450_480_638_690_972_White_et_al_2020_Madagascar_and_surrounding_islands_FCRR.pdf?Signature=AgaeGD2Mv7hctEu%2BYI1IHlhAsQc%3D&Expires=1670592067&AWSAccessKeyId=ASIAWJJOTCLMPJHMKW7W&x-amz-security-token=IQoJb3JpZ2luX2VjEB0aCXVzLXdlc3QtMiJHMEUCIQCNK8cYTCbSyc%2B9AC9sxK44sxNovl9Kh/sgobjv3Jac%2BwIgbh7QlsTWrLPqDQ2eeCRNGo60X%2BhMt4J9zGFxwu6x2H0qzAQIRhADGgw0MzIyNzg4MDMxNjAiDHR5gIcrG23s3%2BJPcSqpBCTYfzrZCOagNotLXzfRhDiwgp7y77nJ8t9TqB1JFQZoK/kyBZI7SgseG1vBA5GBdT8oK5MM6/5nH%2BqrsWGoZqgHTqBdDd8JWGHpbYYjEj600RQMEFoUHouDCompkEccwS4l0jxM97PLBqX%2BB8BsLfykZRuep5Jwp2lKdk8lzONDQTSkpsNgM%2Bj5yHjPPfkOE0sxc68S4Ovcr4maEDyWheLBEdZ9HsmuzuIj5NHV8PX3MVfwJ7A4hyEIhQtCX9uUoJ3yNqMNV/14%2BWEaopxSY4bli04EUmCFIRjp2ctT9qGBGvhDziWDcagn35dS%2BSb/cckPaUSsLI/OOETKcgNszGWqCR%2BznbwrO5uB6DJ4bZNx7u7ZKrszCYIR6DSsUQdxZt0cppIMuIAb5TVk9ALRDMAdwOhHJE2DA%2BkSo9rHnhn2QmoaXIgiyzMuiPs/97Xx5OYQoxbxPz1y2B1mJ8DvF8xFL2nZ9YfHotsCj6dIsR93AaS76Sj3uzrmJofwPMQjNP6Yuv28l53e8Yrilp9%2BkD1D4ACvW9LQMhKZBa1Rft6BflQhesxEzZTzPgHZDB2PfuoMxR7HwRA%2B50jgL8EcAzEKN4WB7Es3k5lTu%2BWxUEsJo2WwTMY6Vd5HEN87ECmeWF%2Bn4Whl4/kPsKz0QSJRIr948V/c78NMotTOqfCz3%2B8ZVZZqF%2BP9X9cb5tF6uGyaJSQAUXroydJXMuIz%2BMJqDKfrCo%2Bp2fKe6F0wtMTHnAY6qQH6u5FifccNsnvtnehEgZOYn2G5NZr8DvVt/fSFpRZfLV5dWybGMMRJosLSy/ZSsV7hEwH76dVZmz9AiVbtngrYj4RMSFOPT23s1%2B8Z9h%2BQDzNhagxY29Kl8/KTshihz6JRJqE/0usT972uaGy30s7HYINnl6qsQBx3lQJ%2BJpvlCcdJFFLdczQBKx3dQseXM6GyDRQtDZBJMNijUWd1N5H2hvPho4gaZwdz


# Load SAU data for Mauritius
sau_mauritius <- read_csv("data/raw-data/sea-around-us/sau-mauritius-12.8.22.csv")

# Filter for most recent year, SSF, and landings only (no discards)
sau_mauritius <- filter(sau_mauritius, year == max(sau_mauritius$year))
sau_mauritius <- filter(sau_mauritius, catch_type == "Landings")
temp <- filter(sau_mauritius, fishing_sector == "Artisanal")
temp2 <- filter(sau_mauritius, fishing_sector == "Subsistence")
sau_mauritius <- bind_rows(temp, temp2)

# Select relevant columns (country, scientific name, and landings)
sau_mauritius <- select(sau_mauritius, area_name, scientific_name, tonnes)

# Aggregate landings by taxon
sau_mauritius <- aggregate(. ~ area_name + scientific_name, data = sau_mauritius, sum)

# Remove unidentified taxa
sau_mauritius <- sau_mauritius[-str_which(sau_mauritius$scientific_name, "identified"), ]

# Remove non-chordata
sau_mauritius <- sau_mauritius[!(sau_mauritius$scientific_name %in% c("Octopodidae", "Octopus cyanea", "Octopus vulgaris", "Panulirus")),]

# Replace landings with ratios
sau_mauritius$ratio <- sau_mauritius$tonnes / sum(sau_mauritius$tonnes)
sau_mauritius <- select(sau_mauritius, !tonnes)

# Get worms results
worms_mauritius <- worms_taxa2(sau_mauritius$scientific_name)

# Combine worms with catch data
mauritius <- bind_cols(worms_mauritius, sau_mauritius)
mauritius <- select(mauritius, valid_name:genus, ratio)




##### Comoros #####




##### * CM SAU #####

## The most recent catch reconstruction for Comoros is:
##    White, R., Noël, S.-L., Christ, H., Relano, V., Sicnawa, F., & Tsui, G. (2020). 
##        Madagascar and smaller islands in the Western Indian Ocean: Updated catch 
##        reconstructions for 2011-2018. In B. Derrick, M. Khalfallah, V. Relano, D. 
##        Zeller, & D. Pauly (Eds.), Updating to 2018 the 1950-2010 marine catch 
##        reconstructions of the Sea Around Us: Part I - Africa, Antarctica, Europe, 
##        and the North Atlantic. (Vol. 5, pp. 81–99). Fisheries Center, University of 
##        British Columbia. https://sau-technical-reports.s3.amazonaws.com/174_175_251_252_450_480_638_690_972_White_et_al_2020_Madagascar_and_surrounding_islands_FCRR.pdf?Signature=AgaeGD2Mv7hctEu%2BYI1IHlhAsQc%3D&Expires=1670592067&AWSAccessKeyId=ASIAWJJOTCLMPJHMKW7W&x-amz-security-token=IQoJb3JpZ2luX2VjEB0aCXVzLXdlc3QtMiJHMEUCIQCNK8cYTCbSyc%2B9AC9sxK44sxNovl9Kh/sgobjv3Jac%2BwIgbh7QlsTWrLPqDQ2eeCRNGo60X%2BhMt4J9zGFxwu6x2H0qzAQIRhADGgw0MzIyNzg4MDMxNjAiDHR5gIcrG23s3%2BJPcSqpBCTYfzrZCOagNotLXzfRhDiwgp7y77nJ8t9TqB1JFQZoK/kyBZI7SgseG1vBA5GBdT8oK5MM6/5nH%2BqrsWGoZqgHTqBdDd8JWGHpbYYjEj600RQMEFoUHouDCompkEccwS4l0jxM97PLBqX%2BB8BsLfykZRuep5Jwp2lKdk8lzONDQTSkpsNgM%2Bj5yHjPPfkOE0sxc68S4Ovcr4maEDyWheLBEdZ9HsmuzuIj5NHV8PX3MVfwJ7A4hyEIhQtCX9uUoJ3yNqMNV/14%2BWEaopxSY4bli04EUmCFIRjp2ctT9qGBGvhDziWDcagn35dS%2BSb/cckPaUSsLI/OOETKcgNszGWqCR%2BznbwrO5uB6DJ4bZNx7u7ZKrszCYIR6DSsUQdxZt0cppIMuIAb5TVk9ALRDMAdwOhHJE2DA%2BkSo9rHnhn2QmoaXIgiyzMuiPs/97Xx5OYQoxbxPz1y2B1mJ8DvF8xFL2nZ9YfHotsCj6dIsR93AaS76Sj3uzrmJofwPMQjNP6Yuv28l53e8Yrilp9%2BkD1D4ACvW9LQMhKZBa1Rft6BflQhesxEzZTzPgHZDB2PfuoMxR7HwRA%2B50jgL8EcAzEKN4WB7Es3k5lTu%2BWxUEsJo2WwTMY6Vd5HEN87ECmeWF%2Bn4Whl4/kPsKz0QSJRIr948V/c78NMotTOqfCz3%2B8ZVZZqF%2BP9X9cb5tF6uGyaJSQAUXroydJXMuIz%2BMJqDKfrCo%2Bp2fKe6F0wtMTHnAY6qQH6u5FifccNsnvtnehEgZOYn2G5NZr8DvVt/fSFpRZfLV5dWybGMMRJosLSy/ZSsV7hEwH76dVZmz9AiVbtngrYj4RMSFOPT23s1%2B8Z9h%2BQDzNhagxY29Kl8/KTshihz6JRJqE/0usT972uaGy30s7HYINnl6qsQBx3lQJ%2BJpvlCcdJFFLdczQBKx3dQseXM6GyDRQtDZBJMNijUWd1N5H2hvPho4gaZwdz

# Load SAU data for Comoros
sau_comoros <- read_csv("data/raw-data/sea-around-us/sau-comoros-12.8.22.csv")

# Filter for most recent year, SSF, and landings only (no discards)
sau_comoros <- filter(sau_comoros, year == max(sau_comoros$year))
sau_comoros <- filter(sau_comoros, catch_type == "Landings")
temp <- filter(sau_comoros, fishing_sector == "Artisanal")
temp2 <- filter(sau_comoros, fishing_sector == "Subsistence")
sau_comoros <- bind_rows(temp, temp2)

# Select relevant columns (country, scientific name, and landings)
sau_comoros <- select(sau_comoros, area_name, scientific_name, tonnes)

# Aggregate landings by taxon
sau_comoros <- aggregate(. ~ area_name + scientific_name, data = sau_comoros, sum)

# Remove unidentified taxa
sau_comoros <- sau_comoros[-str_which(sau_comoros$scientific_name, "identified"), ]

# Replace landings with ratios
sau_comoros$ratio <- sau_comoros$tonnes / sum(sau_comoros$tonnes)
sau_comoros <- select(sau_comoros, !tonnes)

# Get worms results
worms_comoros <- worms_taxa2(sau_comoros$scientific_name)

# Combine worms with catch data
comoros <- bind_cols(worms_comoros, sau_comoros)
comoros <- select(comoros, valid_name:genus, ratio)




##### Save #####

## This section generates nutrient concentrations based on national catch
## taxonomic breakdowns

# Add national catch summaries to a list
wio <- list(comoros = comoros, kenya = kenya, madagascar = madagascar, maldives = maldives, 
  mauritius = mauritius, mayotte = mayotte, mozambique = mozambique,
  reunion = reunion, seychelles = seychelles, tanzania = tanzania)

# Save the list
saveRDS(wio, file = "data/clean-data/04_NationalCatchCompositions.rds")







