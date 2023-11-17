# This script uses Jim Thorson's FishLife package to provide the most accurate estimates available
#   for life history parameters. Unfortunately, it relies on an old version of rfishbase. First,
#   restart R. Then install the old version of fishbase as below. You can then run this code. Once
#   you've run the code, detach rfishbase, update the package back to v. 4 or higher, and restart
#   R again.


# Install old version of fishbase
remotes::install_github( 'ropensci/rfishbase@fb-21.06', force=TRUE )

# NOW RESTART R AND RUN THE CODE FROM THIS POINT

library(rfishbase)
library(FishLife)
library(readr)
library(stringr)

# Import data
species <- read_csv("data/temp-data/02_SpeciesData.csv")

# Empty columns for new life history values
species$Linf_cm <- NA
species$Lmat_cm <- NA
species$Lopt_cm <- NA
species$k <- NA

# Fill in Lopt, Lmat, Linf, and K for all species
for(i in 1:nrow(species)){
  
  # Genus name
  Genus <- str_split_fixed(species$species[i], pattern = " ", n = 2)[1]
  
  # Species name
  Species <- str_split_fixed(species$species[i], pattern = " ", n = 2)[2]
  
  # Predictions of life history parameters for this species (in log space)
  try({ # skip to next iteration of loop if no match
  Predict <- Plot_taxa(Search_species(Genus = Genus, Species = Species)$match_taxonomy, mfrow=c(2,2))
  })
    
  # Predictive mean of Linf for this species
  Linf <- exp(Predict[[1]]$Mean_pred[["Loo"]])
  
  # Predictive mean of Lmat for this species
  Lmat <- exp(Predict[[1]]$Mean_pred[["Lm"]])
  
  # Predictive mean of M for this species
  M <- exp(Predict[[1]]$Mean_pred[["M"]])

  # Predictive mean of K for this species
  K <- exp(Predict[[1]]$Mean_pred[["K"]])
  
  # Calculate Lopt (Beverton 1992)
  Lopt <- Linf * (3 / (3 + M / K))
  
  # Save life history values to species
  species$Linf_cm[i] <- Linf
  species$Lmat_cm[i] <- Lmat
  species$Lopt_cm[i] <- Lopt
  species$k[i] <- K
  
}


# Save species
write.csv(species, file = "data/clean-data/02_SpeciesData.csv",
  row.names = FALSE)

# Detach rfishbase
detach(package:rfishbase, unload = TRUE)

# Update rfishbase to the current version
install.packages("rfishbase")

# NOW RESTART R BEFORE RUNNING ANY OTHER CODE THAT USES FISH BASE


