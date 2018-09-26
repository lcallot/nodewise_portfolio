# Repliation material for: "A Nodewise Regression Approach to Estimating Large Portfolios"
### Mehmet Caner, Esra Ulasan, Laurent Callot, A. Özlem Önder

## Simulation results

To replicate the simulation results in the paper two steps are necessary.

 1. Running the simulations and saving summary statistics. Simulations are launched using the `sim_scripts/sim_setup_script.R` file. Upon completion of the simulation, a `rds` object containing a set of statistics is saved.
 2. Producing the figures in the paper. The figures in the paper are generated using the `sim_scripts/sim_figs.Rnw` Rmarkdown file.


## Empirical Application

The material is located under the `empirical_application` folder and is divided between 4 folders:

  1. `application`: contains a **knitr** file, **application.Rnw**, to be compiled to generate the tables in the paper. **Ndw.Rda** contains all empirical results in the same place.
  2.  `results`: contains **XXX.Rda** files which show all results of different time periods separately.
  3. `code`: contains 4 files. **roll.R** contains the code for Nodewise regression and other high dimensional covariance matrix estimators and the rolling window forecast analysis. **roll_script.R** contains the global settings for the empirical application and runs the experiments. The other two files contain constants and data generation for out-of-sample forecast.
  4. `data`: 5 files called **snp_returns_XXX.rds** contain the asset and risk free rate for 5 different time periods. 5 files called **snp_excess_XXX.rds** contain the excess returns for 5 time periods.
