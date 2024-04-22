#!/bin/env/python

"""
Functions for outputting fit results
"""

import pandas as pd
import itfit
#import bayes
from re import sub

# iterative model fitting
# per guidelines in GAP manual

def params_df(fitlist, drop=["lastfit"]):
    """
    Take a list of fits, return pandas dataframe containing parameters.
    If the list contains original dataArray objects (as opposed to .lastfit),
    the source filename will be in there, too.
    """
    return pd.DataFrame.from_records([itfit.params_get(fit) for fit in fitlist]).drop(drop, axis=1, errors='ignore')

def coords_df(fitdata):
    "return a dataframe of the fitted model coordinates"
    try:
        idx_col  = {val: key for key, val in fitdata.lastfit.columnIndex.items() if val is not None}
        col_ids = [None] * len(idx_col)
        for key, val in idx_col.items():
            col_ids[key] = val
        # extract fitted coords
        coords = pd.DataFrame(fitdata.lastfit).T
        # set colnames
        coords.columns = col_ids
        return coords
    except:
        # for a Bayesian fit, do something else
        coords = pd.DataFrame.from_dict(
            {
                'X': fitdata[0],
                'Y': fitdata.mcmcfit["best_model"]
            }
        )
        return coords

def save_lastfit(fitdata, patt="\.dat", repl="_fitted.dat"):
    """
    Take a fitted dataArray and save the fitted data as a sister file alongside
    the raw data, using the passed regex and replacement as the filename.
    Return the output filename.
    """
    outfile  = sub(patt, repl, fitdata.name)
    with open(outfile, 'w') as handout:
        coords_df(fitdata).to_csv(handout, sep='\t', index=False)