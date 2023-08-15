---
output: html_document
---
# NutrientYieldsWIO

Nutrient yields and fisheries underperformance in the Western Indian Ocean

## Authors

- Bryan P. Galligan S.J. (bgalligan@jesuits.org)

- Timothy R. McClanahan (tmcclanahan@wcs.org)

## Overview

This is the data management and analysis workflow for a study of nutrient yields from coral reef artisanal fisheries in the Western Indian Ocean. It includes 20 years of detailed catch monitoring data with locations in Kenya, Mozambique, and Madagascar, as well as national catch composition reconstructions for 10 WIO exclusive economic zones (EEZs). These data are available for confirmatory purposes only. For any other uses, original data providers must be contacted directly.

## Instructions

### For data access

Processed data files are available in the [`data/clean-data`](https://github.com/bryanpgalligan/NutrientYieldsWIO/tree/main/data/clean-data) subfolder. Files are numbered to indicate which script produces them. See below for metadata and information on data providers.

### To reproduce the analysis

Scripts are designed to be run in numeric order, beginning with `01_ScriptTitle.R`.

## Built with

- R version 4.2.2 (2022-10-31) -- "Innocent and Trusting"

- RStudio 2023.03.0+386 "Cherry Blossom" Release (3c53477afb13ab959aeb5b34df1f10c237b256c3, 2023-03-09) for macOS

### R packages

- corrplot (v0.92)

- dplyr (v1.1.2)

- FishLife (v2.0.1)

- forcats (v1.0.0)

- ggeffects (v1.2.1)

- ggmap (v3.0.2)

- ggplot2 (v3.4.2)

- ggpubr (v0.6.0)

- ggsci (v3.0.0)

- glmmTMB (v1.1.7)

- randomForest (v4.7-1.1)

- readr (v2.1.4)

- readxl (v1.4.2)

- rfishbase (v4.1.1)

### Software citations

Boettiger C, Lang DT, Wainwright PC. (2012). rfishbase: exploring, manipulation, and visualizing FishBase data from R. _Journal of Fish Biology_, _81_(6), 2030--2039. <https://doi.org/10.1111/j.1095-8649.2012.03464.x>

Brooks ME, Kristensen K, van Benthem KJ, Magnusson A, Berg CW, Nielsen A, Skaug HJ, Maechler M, Bolker BM. (2017). glmmTMB balances speed and flexibility among packages for zero-inflated generalized linear mixed modeling. _The R Journal_, _9_(2), 378--400. <https://doi.org/10.32614/RJ-2017-066/>

Kahle D, Wickham H. (2013). ggmap: Spatial visualization with ggplot2. _The R Journal_, _5_(1), 144--161. <https://journal.r-project.org/archive/2013-1/kahle-wickham.pdf>

Kassambara A. (2023). _ggpubr: "ggplot2" based publication ready plots_. R package version 0.6.0. <https://CRAN.R-project.org/package=ggpubr>

Liaw A, Wiener M. (2002) Classification and regression by randomForest. _R News_, _2_(3), 18--22. <https://CRAN.R-project.org/doc/Rnews/>

Lüdecke D. (2018). ggeffects: Tidy data frames of marginal effects from regression models. _Journal of Open Source Software_, _3_(26), 772. <https://doi.org/10.21105/joss.00772>

Posit Team. (2023). _RStudio: Integrated development environment for R_. Posit Software, PBC (Boston, MA). <https://www.posit.co/>

R Core Team. (2022). _R: A language and environment for statistical computing_. R Foundation for Statistical Computing. <https://www.R-project.org/>

Thorson JT. (2020). Predicting recruitment density dependence and intrinsic growth rate for all fishes worldwide using a data-integrated life-history model. _Fish and Fisheries_, _21_(2), 237–-251. <https://doi.org/10.1111/faf.12427>

Thorson JT. (2022). _FishLife: Predict life history parameters for any fish_. R package version 2.0.1. <http://github.com/James-Thorson-NOAA/FishLife>

Thorson JT, Munch SB, Cope JM, Gao J. (2017). Predicting life history parameters for all fishes worldwide. _Ecological Applications_, _27_(8), 2262--2276. <https://doi.org/10.1002/eap.1606>

Wei T, Simko V. (2021). _R package 'corrplot': Visualization of a correlation matrix_. R package version 0.92. <https://github.com/taiyun/corrplot>

Wickham H. (2016). _ggplot2: Elegant graphics for data analysis_. Springer-Verlag (New York). <https://ggplot2.tidyverse.org>

Wickham H. (2023). _forcats: Tools for working with categorical variables (factors)_. R package version 1.0.0. <https://CRAN.R-project.org/package=forcats>

Wickham H, Bryan J. (2023). _readxl: Read excel files_. R package version 1.4.2. <https://CRAN.R-project.org/package=readxl>

Wickham H, François R, Henry L, Müller K, Vaughan D. (2023). _dplyr: A grammar of data manipulation_. R package version 1.1.2. <https://CRAN.R-project.org/package=dplyr>.

Wickham H, Hester J, Bryan J. (2023). _readr: Read rectangular text data_. R package version 2.1.4. <https://CRAN.R-project.org/package=readr>

Xiao N. (2023). _ggsci: Scientific journal and sci-fi themed color palettes for "ggplot2"_. R package version 3.0.0. <https://CRAN.R-project.org/package=ggsci>

## Data files

### Processed (clean) data

### Data sources











