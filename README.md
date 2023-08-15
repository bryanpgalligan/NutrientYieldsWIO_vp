# NutrientYieldsWIO

Nutrient production in tropical fisheries depends on biomass-based management

## Overview

Tropical artisanal fisheries are a crucial source of dietary nutrients for millions of people but are threatened by climate change and overfishing. The need to enhance nutrient production from tropical ecosystems to feed the poor could potentially create a new framework for fisheries science and management. Early recommendations have included targeting small fishes and increasing the species richness of fish catches, which would represent a departure from more traditional biomass-based, sustainable yield, and biodiversity conservation approaches. To test these recommendations, we compared the outcomes of biomass-based management with hypothesized factors influencing nutrient density in artisanal fish catches in large data sets from nearshore fisheries in the Western Indian Ocean. Using machine learning tools and catch reconstruction techniques, we found that enhancing nutrient production depends primarily on achieving biomass-based targets. Therefore, some recommendations designed to enhance the nutrient density of catches are expected to reduce nutrient production, overall yields, and biodiversity.

This repository includes software and summary data associated with our study of nutrient yields from tropical artisanal fisheries in the Western Indian Ocean. All R code is included, but full datsasets are not publicly available. Please contact the developer for assistance gaining access to the original (complete) data.

## Data files

### National Catch Compositions

**Filepath**: [data/04_NationalCatchCompositions.rds](https://github.com/bryanpgalligan/NutrientYieldsWIO_vp/blob/main/data/04_NationalCatchCompositions.rds)

**Description**: The main data file containing catch compositions will appear as a list object when opened in R. It includes 10 data frames containing taxonomic break downs of the artisanal nearshore catches of each WIO fishery jurisdiction. They match the 6 *.csv files listed below, which are the only ones original to this study. The other four are direct copies from the [Sea Around Us (SAU)](https://www.seaaroundus.org/data/#/eez) database and are publicly available on the SAU website.

**Variables**: `9`

| Variable Name | Description |
| ------------- | ----------- |
| valid_name | The scientific name for this taxon |
| rank | Taxonomic rank (family, genus, or species) |
| kingdom | Taxonomic kingdom |
| phylum | Taxonomic phylum |
| class | Taxonomic class |
| order | Taxonomic order |
| family | Taxonomic family |
| genus | Taxonomic genus |
| ratio | Proportion of this jurisdiction's artisanal fish catch taken up by this taxon, expressed as a decimal |

#### Kenya Catch Composition

**Filepath**: [data/KenyaCatchComposition.csv](https://github.com/bryanpgalligan/NutrientYieldsWIO_vp/blob/main/data/KenyaCatchComposition.csv)

**Taxa**: `267`

#### Madagascar Catch Composition

**Filepath**: [data/MadagascarCatchComposition.csv](https://github.com/bryanpgalligan/NutrientYieldsWIO_vp/blob/main/data/MadagascarCatchComposition.csv)

**Taxa*: `197`

#### Mozambique Catch Composition

**Filepath**: [data/MozambiqueCatchComposition.csv](https://github.com/bryanpgalligan/NutrientYieldsWIO_vp/blob/main/data/MozambiqueCatchComposition.csv)

**Taxa*: `133`

#### Reunion Catch Composition

**Filepath**: [data/ReunionCatchComposition.csv](https://github.com/bryanpgalligan/NutrientYieldsWIO_vp/blob/main/data/ReunionCatchComposition.csv)

**Taxa*: `46`

#### Seychelles Catch Composition

**Filepath**: [data/SeychellesCatchComposition.csv](https://github.com/bryanpgalligan/NutrientYieldsWIO_vp/blob/main/data/SeychellesCatchComposition.csv)

**Taxa*: `134`

#### Tanzania Catch Composition

**Filepath**: [data/TanzaniaCatchComposition.csv](https://github.com/bryanpgalligan/NutrientYieldsWIO_vp/blob/main/data/TanzaniaCatchComposition.csv)

**Taxa*: `44`

### National Nutrient Losses

**Filepath**: [data/05_NationalNutrientLosses.csv](https://github.com/bryanpgalligan/NutrientYieldsWIO_vp/blob/main/data/05_NationalNutrientLosses.csv)

**Description**: Summary data of nutrient production and lost production for 10 Western Indian Ocean fishery jurisdictions.

**Observations**: `10`

**Variables**: `39`

| Variable Name | Description |
| ------------- | ----------- |
| eez | Name of fishery jurisdiction (exclusive economic zone) |
| lost.yield_total.tons.day | Lost yield due to overfishing in this jurisdiction (tons per day) |
| lost.yield_mean.kg.km2.day | Mean lost yield at overfished sites (kilograms per square kilometer per day) |
| lost.yield_se.mean | Standard error of the mean for lost yield at overfished sites |
| population | Population estimates for this jurisdiction (World Bank and INSEE) |
| children | Population under age 5 in this jurisdiction (UNDESA, 2022) |
| calcium_pdv.100g | Calcium density in the artisanal catch (percent recommended daily allowance for a child 1--3 years old in 100 g serving) |
| iron_pdv.100g | Iron density in the artisanal catch (percent recommended daily allowance for a child 1--3 years old in 100 g serving) |
| omega3_pdv.100g | Omega-3 density in the artisanal catch (percent allowable intake for a child 1--3 years old in a 100 g serving) |
| selenium_pdv.100g | Selenium density in the artisanal catch (percent recommended daily allowance for a child 1--3 years old in 100 g serving) |
| vitamina_pdv.100g | Vitamin A density in the artisanal catch (percent recommended daily allowance for a child 1--3 years old in 100 g serving) |
| zinc_pdv.100g | Zinc density in the artisanal catch (percent recommended daily allowance for a child 1--3 years old in 100 g serving) |
| lost.calcium_dri.day | Total lost servings of calcium due to overfishing (per day) |
| lost.calcium_dri.child.day | Per capita (child) lost servings of calcium due to overfishing (per day) |
| lost.calcium_mean.dri.km2.day | Mean lost child servings of calcium at overfished sites per square kilometer per day relative to biomass at maximum sustained yield (Bmsy) |
| lost.calcium_se.dri | Standard error of the mean lost servings per square kilometer per day |
| lost.iron_dri.day | Total lost servings of iron due to overfishing (per day) |
| lost.iron_dri.child.day | Per capita (child) lost servings of iron due to overfishing (per day) |
| lost.iron_mean.dri.km2.day | Mean lost child servings of iron at overfished sites per square kilometer per day relative to biomass at maximum sustained yield (Bmsy) |
| lost.iron_se.dri | Standard error of the mean lost servings per square kilometer per day |
| lost.omega3_dri.day | Total lost servings of omega-3 due to overfishing (per day) |
| lost.omega3_dri.child.day | Per capita (child) lost servings of omega-3 due to overfishing (per day) |
| lost.omega3_mean.dri.km2.day | Mean lost child servings of omega-3 at overfished sites per square kilometer per day relative to biomass at maximum sustained yield (Bmsy) |
| lost.omega3_se.dri | Standard error of the mean lost servings per square kilometer per day |
| lost.selenium_dri.day | Total lost servings of selenium due to overfishing (per day) |
| lost.selenium_dri.child.day | Per capita (child) lost servings of selenium due to overfishing (per day) |
| lost.selenium_mean.dri.km2.day | Mean lost child servings of selenium at overfished sites per square kilometer per day relative to biomass at maximum sustained yield (Bmsy) |
| lost.selenium_se.dri | Standard error of the mean lost servings per square kilometer per day |
| lost.vitamina_dri.day | Total lost servings of vitamin A due to overfishing (per day) |
| lost.vitamina_dri.child.day | Per capita (child) lost servings of vitamin A due to overfishing (per day) |
| lost.vitamina_mean.dri.km2.day | Mean lost child servings of vitamin A at overfished sites per square kilometer per day relative to biomass at maximum sustained yield (Bmsy) |
| lost.vitamina_se.dri | Standard error of the mean lost servings per square kilometer per day |
| lost.zinc_dri.day | Total lost servings of zinc due to overfishing (per day) |
| lost.zinc_dri.child.day | Per capita (child) lost servings of zinc due to overfishing (per day) |
| lost.zinc_mean.dri.km2.day | Mean lost child servings of zinc at overfished sites per square kilometer per day relative to biomass at maximum sustained yield (Bmsy) |
| lost.zinc_se.dri | Standard error of the mean lost servings per square kilometer per day |
| nutrient.density_coef | A combined nutrient density coefficient between 0--1, with 1 representing a full serving of all six nutrients in a 100 g serving of fish |
| lost.dri_km2.day | Mean lost nutrient production for all sites in this jurisdiction (recommended intake of six nutrients per square kilometer per day) |
| lost.dri_se | Standard error of the mean (servings per square kilometer per day) |

### Data citations

For a complete list of data sources and citations, please see the supplemental information to the published manuscript. Data sources directly reported and mentioned above, are:

INSEE. (2023). Estimations de population (résultats provisoires arrêtés fin 2022). <https://www.insee.fr/fr/statistiques/fichier/1893198/estim-pop-nreg-sexe-gca-1975-2023.xls>

Pauly D, Zeller D, Palomares MLD. (2020). Sea Around Us concepts, design, and data. <https://seaaroundus.org>

UNDESA. (2022). _World population prospects: The 2022 revision_ [Data portal]. United Nations, Department of Economic and Social Affairs, Population Division. <https://population.un.org/dataportal/>

World Bank. (2023). World Bank open data. <https://data.worldbank.org/>

## Built with

- R version 4.2.2 (2022-10-31) -- "Innocent and Trusting"

- RStudio 2023.06.1+524 "Mountain Hydrangea" Release (547dcf861cac0253a8abb52c135e44e02ba407a1, 2023-07-06) for macOS

### R packages

- corrplot (v0.92)

- dplyr (v1.1.2)

- FishLife (v2.0.1)

- forcats (v1.0.0)

- ggeffects (v1.2.3)

- ggmap (v3.0.2)

- ggplot2 (v3.4.2)

- ggpubr (v0.6.0)

- ggsci (v3.0.0)

- gtable (v0.3.3)

- randomForest (v4.7-1.1)

- readr (v2.1.4)

- readxl (v1.4.3)

- rfishbase (v4.1.2)

- scales (v1.2.1)

- sf (v1.0-14)

- stringdist (v0.9.10)

- stringr (v1.5.0)

- taxize (v0.9.100)

- tidyr (v1.3.0)

- worrms (v0.4.3)

### Software citations

Boettiger C, Lang DT, Wainwright PC. (2012). rfishbase: exploring, manipulation, and visualizing FishBase data from R. _Journal of Fish Biology_, _81_(6), 2030--2039. <https://doi.org/10.1111/j.1095-8649.2012.03464.x>

Chamberlain S, Szocs E. (2013). taxize - taxonomic search and retrieval in R. _F1000Research_, _2_, 191. <https://f1000research.com/articles/2-191/v2>

Chamberlain S, Szoecs E, Foster Z, Arendsee Z, Boettiger C, Ram K, Bartomeus I, Baumgartner J, O'Donnell J, Oksanen J, Greshake Tzovaras B, Marchand P, Tran V, Salmon M, Li G, Grenié M. (2020). _taxize: Taxonomic information from around the web_. R package version 0.9.98. <https://github.com/ropensci/taxize>

Chamberlain S, Vanhoorne B. (2023). _worrms: World Register of Marine Species (WoRMS) Client_. R package version 0.4.3. <https://CRAN.R-project.org/package=worrms>

Kahle D, Wickham H. (2013). ggmap: Spatial visualization with ggplot2. _The R Journal_, _5_(1), 144--161. <https://journal.r-project.org/archive/2013-1/kahle-wickham.pdf>

Kassambara A. (2023). _ggpubr: "ggplot2" based publication ready plots_. R package version 0.6.0. <https://CRAN.R-project.org/package=ggpubr>

Liaw A, Wiener M. (2002) Classification and regression by randomForest. _R News_, _2_(3), 18--22. <https://CRAN.R-project.org/doc/Rnews/>

Lüdecke D. (2018). ggeffects: Tidy data frames of marginal effects from regression models. _Journal of Open Source Software_, _3_(26), 772. <https://doi.org/10.21105/joss.00772>

Pebesma E. (2018). Simple features for R: Standardized support for spatial vector data. _The R Journal_, _10_(1), 439--446. <https://doi.org/10.32614/RJ-2018-009>

Pebesma E, Bivand R. (2023). _Spatial data science: With applications in R_. Chapman and Hall/CRC. <https://doi.org/10.1201/9780429459016>

Posit Team. (2023). _RStudio: Integrated development environment for R_. Posit Software, PBC (Boston, MA). <https://www.posit.co/>

R Core Team. (2022). _R: A language and environment for statistical computing_. R Foundation for Statistical Computing. <https://www.R-project.org/>

Thorson JT. (2020). Predicting recruitment density dependence and intrinsic growth rate for all fishes worldwide using a data-integrated life-history model. _Fish and Fisheries_, _21_(2), 237–-251. <https://doi.org/10.1111/faf.12427>

Thorson JT. (2022). _FishLife: Predict life history parameters for any fish_. R package version 2.0.1. <http://github.com/James-Thorson-NOAA/FishLife>

Thorson JT, Munch SB, Cope JM, Gao J. (2017). Predicting life history parameters for all fishes worldwide. _Ecological Applications_, _27_(8), 2262--2276. <https://doi.org/10.1002/eap.1606>

van der Loo M. (2014). The stringdist package for approximate string matching. _The R Journal_, _6_, 111--122. <https://CRAN.R-project.org/package=stringdist>

Wei T, Simko V. (2021). _R package 'corrplot': Visualization of a correlation matrix_. R package version 0.92. <https://github.com/taiyun/corrplot>

Wickham H. (2016). _ggplot2: Elegant graphics for data analysis_. Springer-Verlag (New York). <https://ggplot2.tidyverse.org>

Wickham H. (2022). _stringr: Simple, consistent wrappers for common string operations_. R package version 1.5.0. <https://CRAN.R-project.org/package=stringr>

Wickham H. (2023). _forcats: Tools for working with categorical variables (factors)_. R package version 1.0.0. <https://CRAN.R-project.org/package=forcats>

Wickham H, Bryan J. (2023). _readxl: Read excel files_. R package version 1.4.2. <https://CRAN.R-project.org/package=readxl>

Wickham H, François R, Henry L, Müller K, Vaughan D. (2023). _dplyr: A grammar of data manipulation_. R package version 1.1.2. <https://CRAN.R-project.org/package=dplyr>

Wickham H, Hester J, Bryan J. (2023). _readr: Read rectangular text data_. R package version 2.1.4. <https://CRAN.R-project.org/package=readr>

Wickham H, Pedersen T. (2023). _gtable: Arrange 'grobs' in tables_. R package version 0.3.3. <https://CRAN.R-project.org/package=gtable>

Wickham H, Seidel D. (2022). _scales: Scale functions for visualization_. R package version 1.2.1. <https://CRAN.R-project.org/package=scales>

Wickham H, Vaughan D, Girlich M. (2023). _tidyr: Tidy Messy Data_. R package version 1.3.0. <https://CRAN.R-project.org/package=tidyr>

Xiao N. (2023). _ggsci: Scientific journal and sci-fi themed color palettes for "ggplot2"_. R package version 3.0.0. <https://CRAN.R-project.org/package=ggsci>