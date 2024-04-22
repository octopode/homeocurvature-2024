#!/bin/env/python

import numpy as np
import statsmodels.api as sm
import scipy
#from rebin import rebin

# components of MCG (Modified Caille Gaussian) model
# default values are from Pabst spreadsheet
def fC(q, sigC=4, rhor=-1):
    "Scalar function for core form factor"
    return np.sqrt(2*np.pi) * sigC * rhor * np.exp(-((sigC*q)**2)/2)
    
def fH(q, zH=18, sigH=3):
    "Scalar or vector function for headgroup form factor"
    return np.sqrt(2*np.pi) * sigH * np.exp(-((sigH*q)**2)/2) * np.cos(q*zH)
    
def sf(q, d, eta1, N):
    """
    Vector function for MCT lamellar structure factor
    Non-integer numbers of stacks are calculated as a linear combination of results for the next lower and higher values.
    (see https://www.sasview.org/docs/user/models/lamellar_hg_stack_caille.html)
    """
    Nlo = int(N)
    Nhi = Nlo+1
    # these are all coefficients inside the sum
    cosomg = np.array([(N-k)*np.cos(k*q*d) for k in range(1, Nhi)])
    expfac = np.array([np.exp(-1 * (d/(2*np.pi))**2 * q**2 * eta1 * np.euler_gamma) for k in range(1, Nhi)])
    pikfac = np.array([(np.pi*k)**(-1 * (d/(2*np.pi))**2 * q**2 * eta1) for k in range(1, Nhi)])
    sumand = cosomg*expfac*pikfac
    # linear combo
    return (Nhi-N)*(Nlo + 2 * sum(sumand[:-1])) + (N-Nlo)*(Nhi + 2 * sum(sumand))
    
def sf_poly(q, d, eta1, N):
    "Polydisperse MCT SF per Frühwirth et al. 2004"
    # N MUST > 2
    # polydispersity
    # number of N-vals sampled
    N = max(N, 2)
    nsf = int(30/np.log(N))
    # truncate to next odd number
    nsf += not nsf%2
    # standard deviation
    if N >= 5: stdev = np.sqrt(N)
    else: stdev = (N-1)/2
    # the array of N-values to sample
    Ns = np.linspace(N-2*stdev, N+2*stdev, nsf)
    weights = np.array([(1/(stdev*np.sqrt(2*np.pi))) * np.exp(-1*(Ni-N)**2/(2*stdev**2)) for Ni in Ns])
    # normalize the weights
    weights = weights/sum(weights)
    strfacs = np.array([sf(q, d, eta1, N=Ni) for Ni in Ns])
    return np.dot(weights, strfacs)
    
def mcg_full(q, N, d, eta1, Nu, sigH, zH, sigC, rhor, scal, offs, vprof=[1], sff=sf_poly):
    "Full lamellar model. np.array vector function for use with jscatter"
    fmfac = fC(q, sigC, rhor) + 2*fH(q, zH, sigH)
    stfac = sff(q, d, eta1, N)
    iq = (1/q**2)*((fmfac**2) * stfac + Nu*(fmfac**2))
    # convolve by PSF before applying scaling
    # this order matters
    # I guess this smearing is good enough for now...
    iq = np.convolve(iq, vprof, mode="same")
    iq = scal * iq + offs
    return iq
    
def mcg_full_georg(q, N, d, eta1, Nu, sigH, zH, sigC, rhor, scal, offs, vprof):
    gamma = np.euler_gamma
    nq = len(q)
    N = int(N)
    
    # FL = sqrt(2*!pi)*par[2]*par[3]*exp(-(par[2]*q)^2/2.)
    FL = np.sqrt(2*np.pi)*sigC*rhor*np.exp(-0.5*(sigC*q)**2)
    # FH = sqrt(2*!pi)*par[1]*exp(-(par[1]*q)^2/2.)*cos(q*par[0])
    FH = np.sqrt(2*np.pi)*sigH*np.exp(-0.5*(sigH*q)**2)*np.cos(q*zH)
    formf = 2*FH + FL
    
    #k = findgen(Nlam-1)+1
    k = np.array(range(N-1)) + 1
    # assuming this is an outer product - JRW
    # it also looks like this is an elementwise cosine
    # has N rows, k cols
    #cosomg = cos(k ## q*d)
    cosomg = np.cos(np.outer(k, q*d))
    #fac = (d/2/np.pi)^2*q^2*eta1
    fac = (2/(2*np.pi))**2 * q**2 * eta1
    exp_fac1 = np.exp(-1*fac*gamma) # this line was trivial to translate
    # pip install rebin
    # from rebin import rebin
    #exp_fac1 = rebin(exp_fac1,nq,N-1)
    exp_fac1 = np.tile(exp_fac1, (N-1,1)) # think this does the job but might need transposing
    #exp_fac2 = transpose(rebin(np.pi*k,N-1,nq))
    exp_fac2 = np.tile(np.pi*k, (nq,1)).transpose() # same shape as fac1 (N x nq)
    # exp_fac3 = exp_fac2^(-rebin(fac,nq,Nlam-1))
    # again, should be elementwise
    exp_fac3 = exp_fac2**(-1*np.tile(fac, (N-1,1))) # N x nq
    prod = cosomg*exp_fac1*exp_fac3 # N x nq
    deltan = N-k
    #s0 = deltan ## prod
    # not sure about this!
    # should end up 1 x nq
    s0 = np.dot(deltan, prod)
    sq = N + 2*s0

    formf_sqrt = formf**2
    I1 = formf_sqrt*sq + Nu*formf_sqrt
    IMCT = I1/q**2
    IMCT = np.convolve(IMCT, np.flip(vprof), mode="same")
    
    return scal * IMCT + offs
    
def sf_only(q, **kwargs):
    "Vector function for lamellar structure factor. Includes Nu term."
    return kwargs["scal"] * np.array([(sf(x, **kwargs) + kwargs["Nu"]/(x**2)) for x in q]) + kwargs["offs"]

# d-guessing functions
def calcd(qs):
    "Helper function takes a list of peak qs, returns mean lamellar d-spacing"
    # in future, better to return list than mean?
    return np.mean([(i+1)*2*np.pi/j for i, j in enumerate(qs)])

def guessd(q, iq, top_n=2):
    """
    Picks up the top_n most intense local maxima.
    Simple method, but sensitive to shoulders.
    """
    peaks = [(q[i], iq[i]) for i in range(len(q))[1:-1] if (iq[i] > iq[i-1]) and (iq[i] > iq[i+1])]
    peaks = sorted(peaks, key=lambda coords: coords[1], reverse=True)[:top_n]
    peaks, _ = list(zip(*peaks))
    return calcd(peaks)

def guessd_lowess(q, iq, frac_smooth=0.01, **kwargs):
    """
    Perform LOWESS smoothing, then detect peaks to guess d-spacing for a lamellar profile.
    Require statsmodels, scipy.
    """
    fitted = sm.nonparametric.lowess(iq, q, frac=frac_smooth, return_sorted=False)
    peaks, _ = scipy.signal.find_peaks(fitted, **kwargs)
    return calcd(peaks)
    
def guessd_cwt(q, iq, top_n=2, widths=np.arange(10, 25), **kwargs):
    """
    Detect peaks to guess d-spacing for a lamellar profile 
    using a contieta1ous wavelet transformation. Require scipy.
    """
    peaks = scipy.signal.find_peaks_cwt(iq, widths, **kwargs)
    # filter to consider only the tallest n peaks
    peaks = sorted(peaks, key=lambda q: iq[q], reverse=True)[:top_n]
    return calcd(peaks)