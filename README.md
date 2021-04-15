# Causal inference for combined RCT and RWD with missing values

This repository contains R notebooks associated with the article (under review): _Treatment effect estimation with missing attributes_.

Three notebooks are available, one for the simulation study and two for the data analysis.

## Simulation study

- [`simulations_combinedNA.Rmd`](https://github.com/imkemayer/combined-incomplete-data/blob/main/simulations_combinedNA.Rmd) reproducing the **simulations and plots** reported in the article.

It needs to be launched either with the code in [`estimators_and_simulations_wo_cw.R`](https://github.com/imkemayer/combined-incomplete-data/blob/main/estimators_and_simulations_wo_cw.R) and [`estimators_wo_cw.R`](https://github.com/imkemayer/combined-incomplete-data/blob/main/estimators_wo_cw.R) that contains the functions (IPSW, CO, AIPSW, and simulation protocol) or with the code in [`estimators_and_simulations.R`](https://github.com/imkemayer/combined-incomplete-data/blob/main/estimators_and_simulations.R) and [`estimators.R`](https://github.com/imkemayer/combined-incomplete-data/blob/main/estimators.R) to also have the Calibration Weighting estimations. Note that the latter requires the function from the [paper](https://arxiv.org/abs/2003.01242) by [Lin Dong](https://lynndung.github.io/about/). Due to an ongoing reviewing process we cannot make the code for the CW estimator available in this repository. You can submit a request to Lin Dong to get access to her implementation of the CW estimator.


## Data analysis

Note that for reproducing the results from this section, an access to the medical data extracted from the [Traumabase registry](http://www.traumabase.eu/en_US) and the [CRASH-2 trial](http://www.crash2.lshtm.ac.uk/) are necessary and these are available upon request from the study investigators.

- [`preprocessCrash2Crash3Traumabase.Rmd`](https://github.com/imkemayer/combined-incomplete-data/blob/main/preprocessCrash2Crash3Traumabase.Rmd) performs the data preprocessing for the joint analysis of CRASH-2 and the Traumabase. It takes as an entry the raw data from each data sets and bind them with proper covariates. The output is the combined data with the raw Traumabase data (with missing values kept). Another similar data frame but with the imputed Traumabase is also produced.

- [`combinedNA_analysis_crash3_traumabase.Rmd`](https://github.com/imkemayer/combined-incomplete-data/blob/main/combinedNA_analysis_crash3_traumabase.Rmd) performs the average treatment effect estimations on the preprocessed data for the separate and joint analyses of CRASH-2 and the Traumabase. The input is the merged table of both the randomized controlled trial and the observational study (corresponding to the output of `preprocessCrash2Crash3Traumabase.Rmd`).

## Acknowledgements

For this work, we have received constructive and fruitful feedback by Bénédicte Colnet, Shu Yang and François-Xavier Ageron.
