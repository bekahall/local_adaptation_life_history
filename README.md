# Local Adaptation of Life History Model

## Mathematica Notebooks

The notebook `analysis.nb` contains all the methods used to perform the analysis of the deterministic island-mainland model of local adaptation. The notebook `plots.nb` is used to generate all the figures in the manuscript, using the output from `analysis.nb` and the Python simulation code described below.

## Python Code

This repo contains the code used to perform the simulations presented in the manuscript. The model implements a Gillespie algorithm to simulate stochastic evolution of a haploid, biallelic population with alleles $A$ and $a$, on an island undergoing seasonal population crashes, with constant migration of $a$ allele individuals from a mainland. A complete description of the model can be found in the manuscript. 

### Requirements

The Python code requires the following packages:

    matplotlib
    pandas
    numpy

### Usage

A single simulation across all parameter values given in `coords.csv` can be perfomed by running the file `main.py` from the terminal:

    main.py run_no

With the input `run_no` setting the run number. `main.py` outputs a CSV file name results_`run_no`.csv in the directory data which gives the gene frequency and total population size before and after the final population crash for each of the parameter values. 

By default, `main.py` sets the initital number of individuals to 1000 for each allele. The default end time is 150 and the parameters are set as $b_A = 0.25$, $b_a = 0.2$, $X_A = 5000$, $X_a = 10000$, $M = 1$. All parameters can be changed within `main.py`.

The role of `main.py` is to run the function `run_simulation` from `simulation.py` across multiple parameter values and output the results in a single file. To run simulations with more customisation, the user may use the function `run_simulation` directly.

To combine the results from multiple simulation runs, the file `calculate_average.py` can be run from the terminal:

    combine_average.py run_no

With the input `run_no` setting the total number of simulation runs that are being combined. 

## CSV Files

`coords.csv` contains all the parameter pairs for $d$ and $\delta$ used to perform the stochastic simulations.

`coarse_coords.csv` contains all the parameter pairs for $d$ and $\delta$ used for the Floquet analysis.

`fine_coords.csv` contains all the parameter pairs for $d$ and $\delta$ used for calculating the local vs. foreign measure of local adaptation and the strength of genetic drift from the deterministic model.

`stoch_data.csv` contains the average gene frequencies and population sizes from the 50 simulation runs, used to create Figure 3 in the manuscript.