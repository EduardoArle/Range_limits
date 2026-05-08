# Range_limits

## Overview

This repository contains the code used to perform ecological niche modelling analyses, SHAP value calculations, GAM analyses, species-level analyses, and figure generation for the manuscript.

## Requirements

Analyses were conducted in R.

Main packages used include:

- xgboost
- mgcv
- shapviz
- sf
- terra
- randomForest
- Boruta
- gratia
- CoordinateCleaner


Analyses were tested on macOS using R 4.3.

No non-standard hardware is required.

Package installation typically requires a few minutes depending on internet connection and system configuration.

Example scripts can typically be run within a few minutes, whereas the complete workflow may require several hours.


## Repository structure

Main workflow scripts are organised sequentially:

- 1–Split_ranges.R
- 2–Select_mammal_species.R
- 3_Taxonomic_harmonisation.R
- 4_Prepare_occurrence_data.R
- 5–Bias_rasters_for_PAs.R
- 6–ENMs_SHAP_values.R
- 7–Variable_correlation.R
- 8–Calculate_metrics.R
- 9–Organise_results_all_species.R
- 10–GAMs.R
- 11–Calculate_slope_results.R
- 12–Boruta_analyses.R


Separate scripts are also provided for reproducing all main and supplementary figures and tables:

- Figure_1.R to Figure_4.R
- Figure_S1.R to Figure_S6.R
- Table_1.R
- Table_S1.R
- Table_S2.R


## Runtime

The complete workflow may require several hours depending on hardware configuration and parallelisation settings.

## Data availability

Occurrence data were obtained from GBIF.
Climate data were obtained from WorldClim.
Species distribution polygons were obtained from the IUCN Red List.

All data sources are cited in the manuscript.

## Code availability

All scripts used for data preparation, ecological niche modelling, SHAP analyses, GAM analyses, species-level analyses, and figure generation are available in this repository.


## Example data

A small example occurrence dataset is provided in:

example_data/Abrothrix_jelskii_occ.csv

This file is included only to demonstrate the workflow structure and script execution. Full analyses rely on the publicly available datasets cited in the manuscript.
