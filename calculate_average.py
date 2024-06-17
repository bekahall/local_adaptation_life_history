import numpy as np
import pandas as pd
import sys
from pathlib import Path

def main():
    """
    Calculates the average gene frequency and population size results across all simulation runs for each set of parameters.

    Outputs file containing these averages titled 'averages.csv' (file_name) in the same director holding original simulation output.

    Requires: 
        - the number of runs (run_no), input when running from the terminal
        - the same file containing all parameter values used for the simulations (coord_file)
        - directory containing simulation output (data_folder). 
          simulation output files should follow the format 'results_#.csv'
    """
    run_no = int(sys.argv[1])

    coord_file = Path("coords.csv")
    data_folder = Path("data/")
    data_folder.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(coord_file)

    df['p0'] = np.zeros(len(df['d']))
    df['N0'] = np.zeros(len(df['d']))
    df['ptau'] = np.zeros(len(df['d']))
    df['Ntau'] = np.zeros(len(df['d']))

    for n in range(1, run_no + 1):
        file_name = Path("data/results_{}.csv".format(n))
        data = pd.read_csv(file_name)
        df['p0'] = df['p0'] + data['p0']
        df['N0'] = df['N0'] + data['N0']
        df['ptau'] = df['ptau'] + data['ptau']
        df['Ntau'] = df['Ntau'] + data['Ntau']

    df['p0'] = df['p0']/run_no
    df['N0'] = df['N0']/run_no
    df['ptau'] = df['ptau']/run_no
    df['Ntau'] = df['Ntau']/run_no

    file_name = 'averages.csv'
    file_location = data_folder / file_name
    df.to_csv(file_location)

main()