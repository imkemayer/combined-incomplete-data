# Causal inference for combined RCT and RWD with missing values

This repository contains R notebooks associated with the article (under review): _Treatment effect estimation with missing attributes_.

Three notebooks are available, one for the simulation study and two for the data analysis.

## Simulation study

- `simulations_combinedNA.Rmd` reproducing the **simulations and plots** reported in the article.

It needs to be launched either with the code in `estimators_and_simulations_wo_cw.R` and `estimators_wo_cw.R` that contains the functions (IPSW, CO, AIPSW, and simulation protocol) or with the code in `estimators_and_simulations.R` and `estimators.R` to also have the Calibration Weighting estimations. Note that the latter requires the function from the [paper](https://arxiv.org/abs/2003.01242) by [Lin Dong](https://lynndung.github.io/about/). Due to an ongoing reviewing process we cannot make the code for the CW estimator available in this repository. You can submit a request to Lin Dong to get access to her implementation of the CW estimator.


## Data analysis
