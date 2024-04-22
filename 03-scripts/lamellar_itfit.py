#!/usr/bin/env python

## This is a "user-friendly" Python script for doing quick MCG fits to SAXS data.
## in addition to installing the dependencies below, you need to make small mod
## to jscatter's dataarray.py. diff my_dataarray.py dataarray.py to see.

## usage: lamellar_itfit.py mysaxs.dat -p > myfit.tsv
## If using -p, you have to close the plot window to execute the next stage of the fit.

## CHANGELOG
## 20230723 JRW: added argparser for command line usage. Functionality is still very basic; fit specs are hardcoded.

import argparse
import sys
import jscatter as jscat
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
import scipy

# components of MCG (Modified Caille Gaussian) model, ported from GAP's IDL
from lamellar import mcg_full, guessd_lowess, guessd_cwt, guessd
from itfit import itfit
import outfit

parser = argparse.ArgumentParser()
parser.add_argument("datfile", type=str, help="path to .dat file")
parser.add_argument("-f", "--fitted", action="store_true", default=False, help="write the fitted data to sister file *_fitted.dat")
parser.add_argument("-p", "--plot-final", action="store_true", default=True, help="plot the fitted data sister to the input file")
parser.add_argument("-i", "--plot-iters", action="store_true", default=False, help="plot each stage of the fit (interactive!)")
parser.add_argument("-q", "--q-domain", type=float, nargs=3, default=[0.05, 0.35, 500], help="q-domain to fit over: min, max, npts")
args = parser.parse_args()

# load experimental data
data = jscat.dA(args.datfile)
# prune (crop and sample) data
data=data.prune(*args.q_domain[:2], int(args.q_domain[2]))

# define fitting iterations using a tuple of dicts of dicts
#fitsteps = (
#    {"fixpar":{'N':25, 'Nu':0.07, 'zH':20, 'eta1':0.05, 'sigH':3, 'sigC':4, 'rhor':-1, 'Nu':0, 'offs':0}, "freepar":{'d':guessd(data[0], data[1]), 'scal':0.0001}, "max_nfev":10},
#    {"fixpar":{'N':25, 'Nu':0.07, 'eta1':0.05, 'sigH':3, 'sigC':4, 'rhor':-1, 'Nu':0, 'offs':0}, "freepar":{'d', 'zH', 'scal'}, "max_nfev":10},
#    {"fixpar":{'N':25, 'Nu':0.07, 'eta1':0.05, 'sigH':3, 'Nu':0, 'offs':0}, "freepar":{'d', 'zH', 'sigC', 'rhor', 'scal'}, "max_nfev":10},
#    {"fixpar":{'N':25, 'eta1':0.05, 'sigH':3, 'Nu':0, 'offs':0}, "freepar":{'d', 'zH', 'sigC', 'rhor', 'Nu', 'scal'}, "max_nfev":10},
#    {"fixpar":{'sigH':3}, "freepar":{'d', 'zH', 'eta1', 'sigC', 'rhor', 'Nu', 'N', 'Nu', 'scal', 'offs'}} # refinement
#)

# assume Caille param = 0
# DOPC chisq ≈ 180
#fitsteps = (
#    {"fixpar":{'N':25, 'Nu':0.07, 'zH':20, 'eta1':0, 'sigH':3, 'sigC':4, 'rhor':-1, 'offs':0.015}, "freepar":{'d':guessd(data[0], data[1]), 'scal':0.01}, "max_nfev":10},
#    {"fixpar":{'N':25, 'Nu':0.07, 'eta1':0, 'sigH':3, 'sigC':4, 'rhor':-1, 'offs':0}, "freepar":{'d', 'zH', 'scal'}, "max_nfev":10},
#    {"fixpar":{'N':25, 'Nu':0.07, 'eta1':0, 'sigH':3, 'offs':0}, "freepar":{'d', 'zH', 'sigC', 'rhor', 'scal'}, "max_nfev":10},
#    {"fixpar":{'N':25, 'eta1':0, 'sigH':3, 'offs':0}, "freepar":{'d', 'zH', 'sigC', 'rhor', 'Nu', 'scal'}, "max_nfev":10},
#    {"fixpar":{'sigH':3, 'eta1':0}, "freepar":{'d', 'zH', 'sigC', 'rhor', 'Nu', 'N', 'Nu', 'scal', 'offs'}} # refinement
#)

# eta1free
#fitsteps = (
#    {"fixpar":{'N':25, 'Nu':0.07, 'zH':20, 'eta1':0, 'sigH':3, 'sigC':4, 'rhor':-1, 'offs':0.015}, "freepar":{'d':guessd(data[0], data[1]), 'scal':0.01}, "max_nfev":10},
#    {"fixpar":{'N':25, 'Nu':0.07, 'eta1':0, 'sigH':3, 'sigC':4, 'rhor':-1, 'offs':0}, "freepar":{'d', 'zH', 'scal'}, "max_nfev":10},
#    {"fixpar":{'N':25, 'Nu':0.07, 'eta1':0, 'sigH':3, 'offs':0}, "freepar":{'d', 'zH', 'sigC', 'rhor', 'scal'}, "max_nfev":10},
#    {"fixpar":{'N':25, 'eta1':0, 'sigH':3, 'offs':0}, "freepar":{'d', 'zH', 'sigC', 'rhor', 'Nu', 'scal'}, "max_nfev":10},
#    {"fixpar":{'sigH':3}, "freepar":{'d', 'zH', 'sigC', 'rhor', 'Nu', 'N', 'Nu', 'scal', 'offs', 'eta1'}} # refinement
#)

# based on DOPC 0 bar result
# **these are the good ones for DOPC!**
#fitsteps = (
#    {"fixpar":{'N':23, 'Nu':0.64, 'eta1':0.086, 'sigH':3, 'sigC':0.615, 'rhor':-10, 'offs':0.018, 'zH':20}, "freepar":{'d':guessd(data[0], data[1]), 'scal':5E-8}, "max_nfev":10},
#    {"fixpar":{'N':23, 'Nu':0.64, 'eta1':0.086, 'sigH':3, 'sigC':0.615, 'rhor':-10, 'offs':0.018}, "freepar":{'d', 'zH', 'scal'}, "max_nfev":10},
#    {"fixpar":{'N':23, 'Nu':0.64, 'eta1':0.086, 'sigH':3, 'offs':0}, "freepar":{'d', 'zH', 'sigC', 'rhor', 'scal'}, "max_nfev":10},
#    {"fixpar":{'N':23, 'eta1':0, 'sigH':3, 'offs':0}, "freepar":{'d', 'zH', 'sigC', 'rhor', 'Nu', 'scal'}, "max_nfev":10},
#    {"fixpar":{'sigH':3}, "freepar":{'d', 'zH', 'sigC', 'rhor', 'Nu', 'N', 'Nu', 'scal', 'offs', 'eta1'}} # refinement
#)

# one more DOPC with proper constraints (rhor was ridiculous)
fitsteps = (
    {"fixpar":{'N':15, 'Nu':2, 'eta1':0.086, 'sigH':3, 'sigC':0.615, 'rhor':-0.7, 'offs':0.018, 'zH':20}, "freepar":{'d':guessd(data[0], data[1]), 'scal':2.7E-7}, "max_nfev":10},
    {"fixpar":{'N':15, 'Nu':2, 'eta1':0.086, 'sigH':3, 'sigC':0.615, 'rhor':-0.7, 'offs':0.018}, "freepar":{'d', 'zH', 'scal'}, "max_nfev":10},
    {"fixpar":{'N':15, 'Nu':2, 'eta1':0.086, 'sigH':3, 'offs':0}, "freepar":{'d', 'zH', 'sigC', 'rhor', 'scal'}, "max_nfev":10},
    {"fixpar":{'N':15, 'eta1':0, 'sigH':3, 'offs':0}, "freepar":{'d', 'zH', 'sigC', 'rhor', 'Nu', 'scal'}, "max_nfev":10},
    {"fixpar":{'sigH':3}, "freepar":{'d', 'zH', 'sigC', 'rhor', 'Nu', 'N', 'Nu', 'scal', 'offs', 'eta1'}} # refinement
)

# for DOPE
## pretty bad but 1500 bar is OK
#fitsteps = (
#    {"fixpar":{'N':23, 'Nu':0.64, 'eta1':0.086, 'sigH':3, 'sigC':0.615, 'rhor':-10, 'offs':0.018, 'zH':20}, "freepar":{'d':guessd(data[0], data[1]), 'scal':1E-8}, "max_nfev":10},
#    {"fixpar":{'N':23, 'Nu':0.64, 'eta1':0.086, 'sigH':3, 'sigC':0.615, 'rhor':-10, 'offs':0.018}, "freepar":{'d', 'zH', 'scal'}, "max_nfev":10},
#    {"fixpar":{'N':23, 'Nu':0.64, 'eta1':0.086, 'sigH':3, 'offs':0}, "freepar":{'d', 'zH', 'sigC', 'rhor', 'scal'}, "max_nfev":10},
#    {"fixpar":{'N':23, 'eta1':0, 'sigH':3, 'offs':0}, "freepar":{'d', 'zH', 'sigC', 'rhor', 'Nu', 'scal'}, "max_nfev":10},
#    {"fixpar":{'sigH':3}, "freepar":{'d', 'zH', 'sigC', 'rhor', 'Nu', 'N', 'Nu', 'scal', 'offs', 'eta1'}} # refinement
#)

# produced 3 in the neighborhood but a bit too thick
#fitsteps = (
#    {"fixpar":{'sigH':3, 'N':30, 'eta1':0.0021, 'offs':0.018, 'Nu':1.56, 'sigC':0.145, 'rhor':-24, 'zH':21.5}, "freepar":{'d':guessd(data[0], data[1]), 'scal':4E-7}, "max_nfev":10},
#    {"fixpar":{'sigH':3, 'N':30, 'eta1':0.0021, 'offs':0.018, 'Nu':1.56, 'sigC':0.145, 'rhor':-24, }, "freepar":{'d', 'zH', 'scal'}, "max_nfev":10},
#    {"fixpar":{'sigH':3, 'N':30, 'eta1':0.0021, 'offs':0.018, 'Nu':1.56}, "freepar":{'d', 'zH', 'sigC', 'rhor', 'scal'}, "max_nfev":10},
#    {"fixpar":{'sigH':3, 'N':30, 'eta1':0.0021, 'offs':0.018}, "freepar":{'d', 'zH', 'sigC', 'rhor', 'Nu', 'scal'}, "max_nfev":10},
#    {"fixpar":{'sigH':3}, "freepar":{'d', 'zH', 'sigC', 'rhor', 'Nu', 'N', 'Nu', 'scal', 'offs', 'eta1'}} # refinement
#)

## Nice and linear, good w constraints, but super thick
#fitsteps = (
#    {"fixpar":{'sigH':3, 'N':30, 'eta1':3.4E-5, 'offs':0.025, 'Nu':3, 'sigC':2.33, 'rhor':-0.68, 'zH':21.5}, "freepar":{'d':guessd(data[0], data[1]), 'scal':4E-7}, "max_nfev":10},
#    {"fixpar":{'sigH':3, 'N':30, 'eta1':3.4E-5, 'offs':0.025, 'Nu':3, 'sigC':2.33, 'rhor':-0.68, }, "freepar":{'d', 'zH', 'scal'}, "max_nfev":10},
#    {"fixpar":{'sigH':3, 'N':30, 'eta1':3.4E-5, 'offs':0.025, 'Nu':3}, "freepar":{'d', 'zH', 'sigC', 'rhor', 'scal'}, "max_nfev":10},
#    {"fixpar":{'sigH':3, 'N':30, 'eta1':3.4E-5, 'offs':0.025}, "freepar":{'d', 'zH', 'sigC', 'rhor', 'Nu', 'scal'}, "max_nfev":10},
#    {"fixpar":{'sigH':3}, "freepar":{'d', 'zH', 'sigC', 'rhor', 'Nu', 'N', 'Nu', 'scal', 'offs', 'eta1'}} # refinement
#)

## Nice and linear, good w constraints, but super thick
#fitsteps = (
#    {"fixpar":{'sigH':3, 'N':30, 'eta1':3.4E-5, 'offs':0.025, 'Nu':3, 'sigC':2.33, 'rhor':-0.68, 'zH':21.5}, "freepar":{'d':guessd(data[0], data[1]), 'scal':4E-7}, "max_nfev":10},
#    {"fixpar":{'sigH':3, 'N':30, 'eta1':3.4E-5, 'offs':0.025, 'Nu':3, 'sigC':2.33, 'rhor':-0.68, }, "freepar":{'d', 'zH', 'scal'}, "max_nfev":10},
#    {"fixpar":{'sigH':3, 'N':30, 'eta1':3.4E-5, 'offs':0.025, 'Nu':3}, "freepar":{'d', 'zH', 'sigC', 'rhor', 'scal'}, "max_nfev":10},
#    {"fixpar":{'sigH':3, 'N':30, 'eta1':3.4E-5, 'offs':0.025}, "freepar":{'d', 'zH', 'sigC', 'rhor', 'Nu', 'scal'}, "max_nfev":10},
#    {"fixpar":{'sigH':3}, "freepar":{'d', 'zH', 'sigC', 'rhor', 'Nu', 'N', 'Nu', 'scal', 'offs', 'eta1'}} # refinement
#)

# too high and too low, in wrong order
# better ones have larger eta1
#fitsteps = (
#    {"fixpar":{'sigH':3, 'N':30, 'eta1':3.4E-5, 'offs':0.025, 'Nu':3, 'sigC':2.33, 'rhor':-0.68, 'zH':20}, "freepar":{'d':guessd(data[0], data[1]), 'scal':4E-7}, "max_nfev":10},
#    {"fixpar":{'sigH':3, 'N':30, 'eta1':3.4E-5, 'offs':0.025, 'Nu':3, 'sigC':2.33, 'rhor':-0.68, }, "freepar":{'d', 'zH', 'scal'}, "max_nfev":10},
#    {"fixpar":{'sigH':3, 'N':30, 'eta1':3.4E-5, 'offs':0.025, 'Nu':3}, "freepar":{'d', 'zH', 'sigC', 'rhor', 'scal'}, "max_nfev":10},
#    {"fixpar":{'sigH':3, 'N':30, 'eta1':3.4E-5, 'offs':0.025}, "freepar":{'d', 'zH', 'sigC', 'rhor', 'Nu', 'scal'}, "max_nfev":10},
#    {"fixpar":{'sigH':3}, "freepar":{'d', 'zH', 'sigC', 'rhor', 'Nu', 'N', 'Nu', 'scal', 'offs', 'eta1'}} # refinement
#)

## flat and too low...but consistent!
#fitsteps = (
#    {"fixpar":{'sigH':3, 'N':30, 'eta1':0.015, 'offs':0.025, 'Nu':3, 'sigC':2.33, 'rhor':-0.68, 'zH':20}, "freepar":{'d':guessd(data[0], data[1]), 'scal':4E-7}, "max_nfev":10},
#    {"fixpar":{'sigH':3, 'N':30, 'eta1':0.015, 'offs':0.025, 'Nu':3, 'sigC':2.33, 'rhor':-0.68, }, "freepar":{'d', 'zH', 'scal'}, "max_nfev":10},
#    {"fixpar":{'sigH':3, 'N':30, 'eta1':0.015, 'offs':0.025, 'Nu':3}, "freepar":{'d', 'zH', 'sigC', 'rhor', 'scal'}, "max_nfev":10},
#    {"fixpar":{'sigH':3, 'N':30, 'eta1':0.015, 'offs':0.025}, "freepar":{'d', 'zH', 'sigC', 'rhor', 'Nu', 'scal'}, "max_nfev":10},
#    {"fixpar":{'sigH':3}, "freepar":{'d', 'zH', 'sigC', 'rhor', 'Nu', 'N', 'Nu', 'scal', 'offs', 'eta1'}} # refinement
#)

# set constant limits
data.setlimit(
    #rhor = [None, 0 ],
    rhor = [-1.5, 0 ], # per GAP manual
    sigC = [-2,  15 ], # ditto, maybe remove?
    eta1 = [0,  1   ], # ditto
    #Nu   = [0,  None],
    #Nu   = [0,  1   ], # per GAP manual, but not sure I implemented it that way?
    # yeah, the way I have Nu set up is different than Georg's Ndiff. His is a scaling
    # factor from 0 (all Bragg) to 1 (all diffuse). Mine is the actual number of
    # uncorrelated bilayers. Technically it should be held < N.
    zH   = [15, 30  ],
    N    = [2,  None],
    offs = [0,  None]
)

# set dynamic constraints
# causes an error in jscatter on execution
# fixed this bug by replacing line 2707 in dataarray.py with:
# nconstrain = sum([not i for i in constrain])
#data.setConstrain(lambda sigC, sigH: (sigC >= sigH))

# a redirect to keep the iterations from coming out on stdout
original = sys.stdout
sys.stdout = sys.stderr
# fit iteratively! (it's a generator function)
for output, plot in itfit(
    data, 
    fitsteps, 
    model  = mcg_full, 
    plot   = True,
    iplot  = args.plot_iters, 
    output = True
): pass
# restore stdout
sys.stdout = original

# save the fitted data if requested
if args.fitted: outfit.save_lastfit(data)

# plot the final fit if requested
if args.plot_final: plot.savefig(args.datfile.replace('.dat', '_fitted.png'))

# write the fitted params to a TSV on stdout
outfit.params_df([data]).to_csv(sys.stdout, sep='\t', index=False)