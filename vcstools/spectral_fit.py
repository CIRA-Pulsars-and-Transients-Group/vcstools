"""
Function used to fit different spectral models to the flux densities of pulsars
"""
import numpy as np
import scipy.optimize as optimization
from iminuit import Minuit
from iminuit.cost import LeastSquares

#matplotlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

import logging
logger = logging.getLogger(__name__)


# Spectral models from Jankowski et al. 2018 and references within
#-----------------------------------------------------------------

def simple_power_law(v, a, b):
    v0 = 1.3e9
    return b*(v/v0)**a

def broken_power_law(v, a1, a2, b, vb):
    v0 = 1.3e9
    x = v / v0
    xb = vb / v0
    y1 = b*x**a1
    y2 = b*x**a2*(xb)**(a1-a2)
    return np.where(x <= xb, y1, y2)

def log_parabolic_spectrum(v, a, b, c):
    v0 = 1.3e9
    x = np.log10( v / v0 )
    return 10**(a*x**2 + b*x + c)

def high_frequency_cut_off_power_law(v, a, b, vc):
    v0 = 1.3e9
    x = v / v0
    xc = vc / v0
    y1 = b*x**(-2) * ( 1 - x / xc )
    y2 = b*x**a
    return np.where(x < xc, y1, y2)

def low_frequency_turn_over_power_law(v, a, b, beta, vc):
    v0 = 1.3e9
    x = v / v0
    xc = v / vc
    return b * x**a * np.exp( a / beta * xc**(-beta) )
    #return b * np.exp( a / beta * xc**(-beta) )

#-----------------------------------------------------------------


def robust_cost_function(f_y, y, sigma_y, k=1.345):
    beta_array = []
    for fi, yi, sigma_i in zip(f_y, y, sigma_y):
        relative_error = (fi - yi)/sigma_i
        if abs(relative_error) < k:
            beta_array.append( 1./2. * relative_error**2 )
        else:
            beta_array.append( k * abs(relative_error) - 1./2. * k**2 )
    return sum(beta_array)

    
def iminuit_fit_spectral_model(freq, flux, flux_err, model=simple_power_law, plot=False, save_name="fit.png"):
    # Model dependent defaults
    if model == simple_power_law:
        # a, b
        start_params = (-1.6, 0.003)
        mod_limits = [(None, 0), (0, None)]
    elif model == broken_power_law:
        # a1, a2, b, vb
        start_params = (-2.6, -2.6, 0.1, 5e8)
        mod_limits = [(None, 0), (None, 0), (0, None), (1e3, 1e9)]
    elif model == log_parabolic_spectrum:
        # a, b, c
        start_params = (-1.6, 1., 1.)
        mod_limits = [(None, 0), (0, None), (0, None)]
    elif model == high_frequency_cut_off_power_law:
        # a, b, vc
        start_params = (-1.6, 1., 1.3e9)
        mod_limits = [(None, 0), (0, None), (1e3, 1e12)]
    elif model == low_frequency_turn_over_power_law:
        # a, b, beta, vc
        start_params = (-2, 1.e1, 1., 2e8)
        mod_limits = [(None, 0), (0, None) , (0., 2.1), (1e3, 1e10)]
    model_str = str(model).split(" ")[1]
    k = len(start_params) # number of model parameters

    # Check if enough inputs
    if len(freq) <= k + 1:
        logger.warn(f"Only {len(freq)} supplied for {model_str} model fit. This is not enough so skipping")
        return 1e9, None, None, None

    # Fit model
    least_squares = LeastSquares(freq, flux, flux_err, model)
    least_squares.loss = "soft_l1"
    m = Minuit(least_squares, *start_params)
    m.limits = mod_limits
    m.scan(ncall=50)
    m.migrad()  # finds minimum of least_squares function
    m.hesse()   # accurately computes uncertainties
    logger.debug(m)


    # display legend with some fit info
    fit_info = [model_str]
    for p, v, e in zip(m.parameters, m.values, m.errors):
        fit_info.append(f"{p} = ${v:.5f} \\pm {e:.5}$")
    fit_info = "\n".join(fit_info)

    # Calculate AIC
    beta = robust_cost_function(model(np.array(freq), *m.values), flux, flux_err)
    aic = 2*beta * 2*k + (2*k*(k+1)) / (len(freq) - k -1)

    if plot:
        fitted_freq = np.linspace(min(freq), max(freq), 100)
        fitted_flux = model(fitted_freq, *m.values)#, v0=np.mean(freq))
        #fitted_flux = model(fitted_freq, *(-1.5, 0.1, 3., 4000000))
        fig, ax = plt.subplots()
        plt.errorbar(np.array(freq) / 1e9, flux, yerr=flux_err, fmt='o', label="Input data", color="orange")
        plt.plot(fitted_freq / 1e9, fitted_flux, 'k--', label=fit_info) # Modelled line
        plt.xscale('log')
        plt.yscale('log')
        ax.get_xaxis().set_major_formatter(ScalarFormatter())
        ax.get_yaxis().set_major_formatter(ScalarFormatter())
        plt.xlabel('Frequency (GHz)')
        plt.ylabel('Flux (Jy)')
        plt.legend()
        plt.savefig(save_name)
        plt.clf()
    return aic, m.parameters, m.values, m.errors


def find_best_spectral_fit(pulsar, freq_all, flux_all, flux_err_all, plot=False):
    # loop over models
    models = [
            [simple_power_law, "simple_power_law"],
            [broken_power_law, "broken_power_law"],
            [log_parabolic_spectrum, "log_parabolic_spectrum"],
            [high_frequency_cut_off_power_law, "high_frequency_cut_off_power_law"],
            [low_frequency_turn_over_power_law, "low_frequency_turn_over_power_law"],
            ]
    aics = []
    fit_results = []
    for model, label in models:
        #curve_fit_spectral_model(freq_all, flux_all, flux_err_all, model=model, plot=True, save_name=f"{pulsar}_{label}_fit.png")
        aic, parameters, values, errors = iminuit_fit_spectral_model(freq_all, flux_all, flux_err_all, model=model, plot=plot, save_name=f"{pulsar}_{label}_fit.png")
        logger.debug(f"{label} model fit gave AIC {aic}.")
        aics.append(aic)
        fit_results.append([parameters, values, errors])
    aici = aics.index(min(aics))
    logger.info(f"Best model for {pulsar} is {models[aici][1]}")
    return models[aici], fit_results[aici]