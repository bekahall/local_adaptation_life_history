import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
from simulation import *

def main():
    """
    Runs a single simulation for each set of parameters given in the chose coordinate file (coord_file). 

    Outputs a csv file named 'results_#.csv' (file_name) into the directory 'data', 
    where # is the run number input when running from the terminal (run_no).
    """
    # Parameters
    nA = 1000
    na = 1000
    
    bA = 0.25
    ba = 0.2
    XA = 5000
    Xa = 10000
    M = 1

    end_time = 150

    run_no = int(sys.argv[1])

    coord_file = Path("coords.csv")
    data_folder = Path("data/")
    data_folder.mkdir(parents=True, exist_ok=True)

    coords = pd.read_csv(coord_file)

    coords['p0'] = np.zeros(len(coords['d']))
    coords['N0'] = np.zeros(len(coords['d']))
    coords['ptau'] = np.zeros(len(coords['d']))
    coords['Ntau'] = np.zeros(len(coords['d']))

    for ind in coords.index:
        delta = coords['delta'][ind]
        d = coords['d'][ind]
        pars = [bA, ba, XA, Xa, d, M]

        [tracknA, trackna, trackt, trackp0, trackN0, trackptau, trackNtau] = run_simulation(nA,na,delta,pars,end_time)
        coords['p0'][ind] = trackp0[-1]
        coords['N0'][ind] = trackN0[-1]
        coords['ptau'][ind] = trackptau[-1]
        coords['Ntau'][ind] = trackNtau[-1]

    file_name = 'results_{}.csv'.format(run_no)
    file_location = data_folder / file_name
    coords.to_csv(file_location)

main()