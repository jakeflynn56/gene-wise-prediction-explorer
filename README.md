# Gene-wise Prediction Explorer
This interactive tool allows you to explore gene-wise prediction using a dataset of your choice. Select the dataset, predictor, and the gene you want to investigate. The app will provide you with visualizations of predictor score distributions, helping you gain an understanding of predictor performance at the gene level. The density plot is estimated with scores from all SNVs.

## Features
- **Dataset Selection:** Users can choose datasets for analysis.
- **Predictor Analysis:** Incorporates two predictors: REVEL and BayesDel.
- **Gene Filtering:** Allows filtering by specific genes to analyze discordant, concordant, and indeterminate intervals.
- **Interactive Visualizations:** Provides density plots and histograms for analysis.
- **Data Export:** Options to download interactive plots as PNG and filtered data as CSV.

## Installation and Usage
- The application can be accessed on [fowler-web](https://fowler-shiny.gs.washington.edu/shiny/fowler-shiny/).
- Alternatively, the application can be ran locally:
    - **Clone the repository:** ```git clone https://github.com/jakeflynn56/gene-wise-prediction-explorer```
    - **Set up the environment:** Ensure you have R and RStudio installed. Dependencies are listed in requirements.R.
    - **Download the required datasets:** See Data section below.
    - **Load the app:** Open app.R in RStudio and run the app.
    - **Explore:** Select your dataset, predictors, and genes to view the visualizations.

## Files in the Repository
- **app.R:** The main R Shiny application script.
- **GenomicDataProcessingScript.sql:** SQL script for creating the SQLite database used for generating density plots. 
- **SNVsConcatenationScript.R:** R script for concatenating SNVs and their BayesDel scores.
- **requirements.R:** Lists all required R packages.
- **README.md:** This file, containing detailed information about the app.

## Repository Structure
```
gene-wise-prediction-explorer/
├── app.R
├── GenomicDataProcessingScript.sql
├── SNVsConcatenationScript.R
├── requirements.R
├── README.md
└── LICENSE
```

## Data
The datasets required for running this application can be found on [Google Drive](https://drive.google.com/drive/folders/1tbs8NvXBmrcviAgPv6Bmt7SWYU3laTTj?usp=sharing). The datasets required for creating the SQLite database can be found [here](https://drive.google.com/drive/folders/1K4LI6ZSsUGBhHoChUtegC8bgCt7hbQlA) for BayesDel and [here](https://drive.google.com/drive/folders/1K4LI6ZSsUGBhHoChUtegC8bgCt7hbQlA) for REVEL. **It is important to note that as of 12/14/2023 the ClinGen SVI Calibration Dataset and ClinVar 2023 Dataset need to be analyzed before the filters will work successfully within the R Shiny application.**

## Authors
Jake Flynn - Initial Work

## License
This project is licensed under the MIT License - see the LICENSE file for details.
