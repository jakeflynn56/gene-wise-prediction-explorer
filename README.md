# Gene-wise Prediction Explorer
This interactive tool allows you to explore gene-wise prediction using a dataset of your choice. Select the dataset, predictor, and the gene you want to investigate. The app will provide you with visualizations of predictor score distributions, helping you gain a deeper understanding of predictor performance at the gene level. The density plot is estimated with scores from all SNVs.

## About
The Gene-wise Prediction Explorer is an interactive R Shiny application designed to facilitate gene-wise prediction analysis using custom datasets. It offers insightful visualizations of predictor score distributions, aiding in the understanding of predictor performance at the gene level.

## Features
- **Dataset Selection:** Users can choose datasets for analysis.
- **Predictor Analysis:** Incorporates two predictors: REVEL and BayesDel.
- **Gene Filtering:** Allows filtering by specific genes to analyze discordant, concordant, and indeterminate intervals.
- **Interactive Visualizations:** Provides density plots for in-depth analysis.
- **Data Export:** Options to download interactive plots as PNG and filtered data as CSV.

## Installation and Usage
- The application can be accessed on [fowler-web](https://fowler-shiny.gs.washington.edu/shiny/fowler-shiny/).
- Alternatively, the application can be ran locally
    - **Clone the repository:** ```git clone [repository URL]```
    - **Set up the environment:** Ensure you have R and RStudio installed. Dependencies are listed in requirements.R.
    - **Load the App:** Open app.R in RStudio and run the app.
    - **Explore:** Select your dataset, predictors, and genes to view the visualizations.

## Files in the Repository
- **app.R:** The main R Shiny application script.
- **GenomicDataProcessingScript.sql:** SQL script for creating the SQLite database used for generating density plots. 
- **SNVsConcatenationScript.R:** R script for concatenating SNVs and their BayesDel scores.
- **requirements.R:** Lists all required R packages.
- **README.md:** This file, containing detailed information about the app.

## Repository Structure
gene-wise-prediction-explorer/
│
├── app.R
├── GenomicDataProcessingScript.sql
├── SNVsConcatenationScript.R
├── requirements.R
├── README.md
│
├── data/
│   ├── sample_data.csv
│   └── bd_no_training_var.csv
│   └── bd_one_star_missense.csv
│   └── revel_no_training_var.csv
│   └── revel_one_star_missense.csv
│   └── clinvar_2023_annoavar_37_parsed_one_star_missense_no_vus_AF_filtered_w_geneid_cleanedup_nogenedup.csv
│   └── clinvar23_no_training_var_AF_fixed_w_geneid_cleaned_up.csv
│   └── Clinvar_2019_dataset_080823.xlsx
│   └── Clinvar_2020_dataset_080823.xlsx
│   └── clinvar_2023_annoavar_37_parsed_one_star_missense_no_vus_AF_filtered.csv
│
└── LICENSE

## Authors
Jake Flynn - Initial Work

## License
This project is licensed under the MIT License - see the LICENSE file for details.
