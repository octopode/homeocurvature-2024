#!/bin/env/python

import matplotlib.pyplot as plt
import lamellar # needed for guess

# iterative model fitting
# per guidelines in GAP manual

def pnames_get(fit, excl=("cov", "func_code", "func_name")):
    "Get list of free param names from jscatter dataArray"
    return [attr for attr in fit.attr if attr not in excl]

def params_get(fit):
    "Extract dict of fitted params from jscatter dataArray"
    return {attr: getattr(fit, attr) for attr in pnames_get(fit)}
    
def gen_fititers(params, fixfree):
    """
    From a dict of parameter values and a list of fixed/free specs, 
    (from configfile) return a list of fititers for itfit.
    """
    # sub the values into the input structure
    first = True
    for iter in fixfree:
        iter["fixpar"] = {par: params[par] for par in iter["fixpar"]}
        if first: 
            iter["freepar"] = {par: params[par] for par in iter["freepar"]}
            first = False
    return fixfree
        
def itfit(data, fititers, iplot=False, plot=False, giveup=False, **kwargs):
    """
    Iterative fitting routine for an efficient search of high-dimensional
    solution space.
    
    data is a jscatter dataArray
    
    ***
    giveup is a switch I added to dataArray.fit() to allow non-converging
    iters. See dataArray.py line 3377.
    ***
    
    fititers has a very specific format:
    
    This is now a generator function that yields the fitted dataArray following
    every fit iter (but not every func evaluation).

    iplot switches interactive plotting
    plot has a plot returned from the generator
    """
    for i, iter in enumerate(fititers):
        # intial run
        if not i:
            data.fit(mapNames={'q':'X'}, giveup=giveup, **iter, **kwargs)
        # subsequent iters
        else:
            try:
                # extract freepar from the previous run
                # .lastfit may be unnecessary
                params_last = params_get(data.lastfit)
            except:
                params_last = params_get(data)
                pass
            # insert those starting values
            iter["freepar"] = {p:params_last[p] for p in iter["freepar"]}
            # fit again
            data.fit(mapNames={'q':'X'}, giveup=giveup, **iter, **kwargs)
        # plot the data and fit
        if plot or iplot:
            plt.cla() # clear the plot axes
            plt.plot(data[0], data[1], color="C0")
            plt.plot(data.lastfit[0], data.lastfit[1], color="C1")
            plt.yscale("log")
            plt.xlabel("q (1/Å)")
            plt.ylabel("I(q) (arb.)")
            plt.title("Fit stage {}".format(i))
        if iplot: plt.show()
        if plot: yield data, plt
        else: yield data