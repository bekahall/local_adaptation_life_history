import matplotlib.pyplot as plt
import numpy as np
import math

def run_simulation(nA,na,delta,pars,end_time,tau=6.28,display=False):
    """
    Runs a single simulation.

    Parameters:
        nA (int): initial number of individuals carrying the 'A' allele
        na (int): initial number of individuals carrying the 'a' allele
        delta (float): proportion of the population that is removed at the end of every season
        pars (list): list of parameters in the following order
                        - bA (float): birth rate of the A allele
                        - ba (float): birth rate of the a allele
                        - XA (int): growth limit of the A allele
                        - Xa (int): growth limit of the a allele
                        - d (float): death rate
                        - M (int, float): migration rate of a alleles onto the island
        end_time (int): number of time steps to run the simulation
        tau (float, default = 6.28): length of a single single
        display (bool, default = False): set whether a plot of the simulation results will be displayed at the end of the simulation

    Returns:
        tracknA (list): list of the number of A allele individuals at each update
        trackna (list): list of the number of a allele individuals at each update
        trackt (list): list of the amount of time passed at each update
        trackp0 (list): list of the gene frequency at the beginning of each season
        trackN0 (list): list of the total population size at the beginning of each season
        trackptau (list): list of the gene frequency at the end of each season
        trackNtau (list): list of the total population size at the end of each season
    """

    # preallocate lists for updating population sizes and times
    dt_max = preallocate(pars)
    length = math.ceil(end_time/dt_max)
    tracknA = np.zeros(length)
    trackna = np.zeros(length)
    trackt = np.zeros(length)
    tracknA[0] = nA
    trackna[0] = na

    i = 1   # index for population sizes and times

    # preallocate lists for updating gene frequencies and total population sizes at the beginning and end of a season
    length = math.ceil(end_time/tau) + 1
    trackp0 = np.zeros(length)
    trackN0 = np.zeros(length)
    trackptau = np.zeros(length)
    trackNtau = np.zeros(length)
    trackp0[0] = nA/(nA+na)
    trackN0[0] = nA + na

    j = 1   # index for gene frequencies and total population sizes

    t = 0
    counter = 0   # will track time since the beginning of the season

    # begin simulation
    while t < end_time:
        # calculate event rates
        bA = calc_bA(nA,na,pars)
        ba = calc_ba(nA,na,pars)
        dA = calc_dA(nA,na,pars)
        da = calc_da(nA,na,pars)
        M = pars[5]

        dt = 1/(bA + ba + dA + da + M)

        # randomly choose an event and update population sizes
        rand = np.random.uniform()

        if rand < bA * dt:
            nA = nA + 1   # birth of A allele individual
        elif rand < (bA + ba) * dt:
            na = na + 1   # birth of a allele individual
        elif rand < (bA + ba + dA) * dt:
            nA = nA - 1   # death of A allele individual
        elif rand < (bA + ba + dA + da) * dt:
            na = na - 1   # death of a allele individual
        elif rand < (bA + ba + dA + da + M) * dt:
            na = na + 1   # migration of a allele individual

        t = t + dt
        counter = counter + dt

        if counter > tau:    # check if the length of a full season has passed
            # save gene frequency and total population size before the crash
            counter = 0
            nT = nA + na
            if nT > 0:
                p = nA/nT
            else:
                p = 0
            trackptau[j] = p
            trackNtau[j] = nT

            # randomly remove a proportion of the population
            deaths = math.floor(delta*nT)
            nAdeaths = np.random.binomial(deaths,p)
            if nAdeaths > nA:
                nAdeaths = nA
                nadeaths = deaths - nA
            elif deaths - nAdeaths > na:
                nadeaths = na
                nAdeaths = deaths - na
            else:
                nadeaths = deaths - nAdeaths
            nA = nA - nAdeaths
            na = na - nadeaths

            # save gene frequency and total population size after the crash
            nT = nA + na
            if nT > 0:
                p = nA/nT
            else:
                p = 0

            trackp0[j] = p
            trackN0[j] = nT
            j = j + 1

        # updation population sizes and times
        tracknA[i] = nA
        trackna[i] = na
        trackt[i] = t
        i = i + 1

    # remove trailing zeroes 
    tracknA = tracknA[:i]
    trackna = trackna[:i]
    trackt = trackt[:i]
    trackp0 = trackp0[:j]
    trackN0 = trackN0[:j]
    trackptau = trackptau[:j]
    trackNtau = trackNtau[:j]

    if display:

        Nlist = tracknA + trackna
        plist = tracknA/Nlist

        fig, ax1 = plt.subplots()

        color = 'darkslateblue'
        ax1.set_xlabel('time',size=16)
        ax1.set_ylabel('population size', color=color,size=16)
        ax1.plot(trackt, Nlist, color=color)
        ax1.tick_params(axis='y', labelcolor=color, labelsize=12)

        ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

        color = 'chocolate'
        ax2.set_ylabel('allele frequency', color=color, size=16)  # we already handled the x-label with ax1
        ax2.plot(trackt, plist, color=color)
        ax2.tick_params(axis='y', labelcolor=color, labelsize=12)

        fig.tight_layout()  # otherwise the right y-label is slightly clipped

        plt.show()

        print(nA,na)
        print(nA + na)
        if (nA + na) > 0:
            print(nA/(nA + na))
    
    return tracknA, trackna, trackt, trackp0, trackN0, trackptau, trackNtau

def calc_bA(nA,na,pars):

    return pars[0] * nA

def calc_ba(nA,na,pars):

    return pars[1] * na

def calc_dA(nA,na,pars):

    return (pars[4] + pars[0] * (nA + na)/pars[2]) * nA

def calc_da(nA,na,pars):

    return (pars[4] + pars[1] * (nA + na)/pars[3]) * na

def preallocate(pars):

    bA = pars[0]
    ba = pars[1]
    XA = pars[2]
    Xa = pars[3]
    d = pars[4]
    M = pars[5]

    X = max(Xa,XA)

    Arate = calc_bA(X,0,pars) + calc_ba(X,0,pars) + calc_dA(X,0,pars) + calc_da(X,0,pars) + M
    arate = calc_bA(0,X,pars) + calc_ba(0,X,pars) + calc_dA(0,X,pars) + calc_da(0,X,pars) + M

    rate = max(Arate,arate)
    return 1/rate