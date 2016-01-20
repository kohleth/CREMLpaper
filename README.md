# The datasets and scripts that accompany the C-REML paper.
## To install C-REML on R.
A developmental version of C-REML is on another github repository and can be assessed by the command ```devtools::install_github("kohleth/spcreml")``` in R.

## Script for the simulation study in the paper.
The folder simStudy contains the script simstudy.R which runs the simulation study. It compares geoR, asreml, and C-REML. This means you need to have all three packages installed on R. Note that asreml is a proprietary package. 

Depending on the processing power of your computer, you might want to adjust the parameter ```Nlist``` and ```Nsim``` accordingly. ```Nlist``` controls the sample size of the generated dataset. ```Nsim``` controls the number of replication in the study.

## Script for the analysis of pH data.
The folder pH_analysis contains the datasets and script for the analysis of the pH dataset.
The script ```analysis.R``` is the main script which contains the analysis code. If you run this on a normal computer with <20 cores, it will take a long time to fit the model. Therefore we have provided the fitted model ```fit1.rda``` so you can just load it onto R ```load("fit1.rda")``` when necessary. 

The prediction grid in this script is done on a much lower resolution than the 50X50m said in the paper, this is so that the file size of the covariate map (```covar.tif```) can be greatly reduced (from >1 GB to <1 MB). 
