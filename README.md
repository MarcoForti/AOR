# README

This repository contains MATLAB code for performing Fuzzy K-Means clustering with L0 penalisation in the membershi matrix. 
The code is implemented to determine the optimal number of clusters and penalty parameter based on the Xie & Beni (XB) index.

## Contents
- `FKM_L0_Main.m`: Script to set-up parameters, options and to execute the functions.
- `FKM.m`: Function to perform the Fuzzy K-Means clustering.
- `FKM_L0.m`: Function to implement FKM-L0 (Non-Exhaustive) clustering with penalty parameter.
- `FKM_L0_Lambda.m`: Function to perform FKM-L0 (Exhaustive) clustering with penalty parameter.
- `XB.m`: Function to compute the Xie & Beni (XB) index for evaluating cluster validity.
- `README.md`: This file explaining the repository contents.
- `Citation.txt`: This file explaining the article to cite.

## Usage
1. **Data Generation/Loading**: 
   - The data must be generated or uploaded.

2. **Set Parameters**:
   - Modify parameters such as `m` (Fuzzyfing parameter), `Max_C` (maximum number of clusters to consider), `Max_iter` (maximum number of iterations), `stand` (standardization option), and `alpha` (parameter for FKM-L0).

3. **Run `FKM_L0_Main.m`**:
   - Execute `FKM_L0_Main.m` to determine the optimal number of clusters (`C`) and penalty parameter (`lambda`) based on the Xie & Beni index.

4. **Interpret Results**:
   - The script will compute the Xie & Beni index for different cluster sizes (`C`) and penalty parameter (`lambda`) values.
   - The optimal number of clusters (`C`) and penalty parameter (`lambda`) that minimize the index are determined.
   - The final Fuzzy K-Means clustering with the optimal parameters is performed, and cluster assignments are obtained.

## Requirements
- MATLAB (R2019a or later)
- Statistics and Machine Learning Toolbox

## Notes
- The code employs parallel computing (`parfor`) for efficiency in evaluating the Xie & Beni index for multiple random starts and cluster sizes.
- Ensure to adjust the data generation part (`X` matrix) or load your own dataset for clustering tasks.

Please feel free to reach out for any questions or improvements to the code. Contributions are welcome!
