# Maize_Prediction_Project
This is the Github repository with the input data and example scripts for the following paper.

## File Overview
- `README.md` - Project overview and file descriptions.
- `Data_Processing.ipynb` - Notebook for preparing and cleaning data used in downstream analyses.
- `Feature_Investigation.ipynb` - Notebook exploring and evaluating candidate input features.
- `Expression_Data_SVR_Prediction.ipynb` - Notebook for support vector regression (SVR) modeling using expression data. 
    - Note: all SVR models were trained and tested the same way, so all that changes between each notebook is the input dataset
- `Label_Investigation.ipynb` - Notebook for inspecting and summarizing prediction labels/targets.
- `Phenotype_Group_Stacked_Bar.ipynb` - Notebook generating stacked bar visualizations for phenotype groups.
- `rrblup_exp.R` - R script for rrBLUP-based analysis on expression-related data.
    - Note: all rrBLUP models were trained and tested the same way, so all that changes between each notebook is the input dataset
