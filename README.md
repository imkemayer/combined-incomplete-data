# Causal inference for combined RCT and RWD with missing values

This repository contains R notebooks associated with the article (under review): _Treatment effect estimation with missing attributes_.

Three notebooks are available, one for the simulation study and two for the data analysis.

## Simulation study

- `simulations_combinedNA.Rmd` reproducing the **simulations and plots** reported in the article.

It needs to be launched with the `estimators_and_simulations-wo-cw.R` code that contains the functions (IPSW, CO, AIPSW, and simulation protocol).

This notebook does not contain the calibration weighting (CW) method as the function used in the paper is the one from [Lin Dong](https://lynndung.github.io/about/) as an implementations of her [research work](https://arxiv.org/abs/2003.01242). Due to an ongoing reviewing process we only publish the results of the simulations and results without the CW method.


## Data analysis
