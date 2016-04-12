import matplotlib.pyplot as plt
from numpy import log, cos, exp, arange
import numpy as np
import gasabsr98
import pandas as pd

# ================================================================ 
# WE ARE FINISHED LOADING AND PRE-PROCESSING OUR ATMOSPHERIC PROFILE.
# EVERYTHING BELOW THIS LINE PERTAINS TO THE FREQUENCY-DEPENDENT
# RADIATIVE TRANSFER CALCULATIONS
#=================================================================

def get_weighting_function_upwelling(f, layers, mu=1.):

    # Get mass extinction for air
    layers['ke_air'] = layers.apply(
        lambda s: gasabsr98.gasabsr98(f,s.Tbar, s.rhowv, s.pbar)[0],
        axis=1
    )

    # Get mass extinction for water vapor
    layers['ke_wv'] = layers.apply(
        lambda s: gasabsr98.gasabsr98(f,s.Tbar, s.rhowv, s.pbar)[1],
        axis=1
    )

    # Get optical depth from extinction
    layers['tau'] = (layers.ke_air*(layers.mass - layers.wvmass)
        + layers.ke_wv*layers.wvmass)

    # Set index to be layer number
    layers = layers.reset_index(drop=True)

    # Compute transmission from top of atmosphere to each layer
    layers['tx_from_top'] = (layers
                    # Take tau(z)
                     .sort_values('Zbar', ascending=False).tau
                    # Compute total optical depths from toa to z for each z
                     .cumsum()
                    # Compute tx = e^(-tau/mu)
                    .apply(lambda tau: np.exp(-tau/mu)))

    # Compute the contribution from each layer to weighting function
    layers['layer_weight'] = -(layers.sort_values('Zbar', ascending=False).tx_from_top
        # Take derivative of transmission
        .diff())

    # Get weight for each layer in m^-1
    layers['weighting_function'] = (layers.layer_weight 
        # by dividing by height in meters
        / layers.sort_values('Zbar', ascending=True).Zbar.diff() * 1000)
    return layers

def get_weighting_function_downwelling(f, layers, mu=1):

    # Get mass extinction for air
    layers['ke_air'] = layers.apply(
        lambda s: gasabsr98.gasabsr98(f,s.Tbar, s.rhowv, s.pbar)[0],
        axis=1
    )

    # Get mass extinction for water vapor
    layers['ke_wv'] = layers.apply(
        lambda s: gasabsr98.gasabsr98(f,s.Tbar, s.rhowv, s.pbar)[1],
        axis=1
    )

    # Get optical depth from extinction
    layers['tau'] = (layers.ke_air*(layers.mass - layers.wvmass)
        + layers.ke_wv*layers.wvmass)

    # Set index to be layer number
    layers = layers.reset_index(drop=True)

    # Compute transmission from top of atmosphere to each layer
    layers['tx_from_bottom'] = (layers
                    # Take tau(z)
                     .sort_values('Zbar', ascending=True).tau
                    # Compute total optical depths from toa to z for each z
                     .cumsum()
                    # Compute tx = e^(-tau/mu)
                    .apply(lambda tau: np.exp(-tau/mu)))

    # Compute the contribution from each layer
    layers['layer_weight'] = (layers.sort_values('Zbar', ascending=False).tx_from_bottom
        # Take derivative of transmission
        .diff())

    # Get weight for each layer in m^-1
    layers['weighting_function'] = (layers.layer_weight
        # by dividing by height in meters
        / -layers.sort_values('Zbar', ascending=False).Zbar.diff() * 1000)

    return layers



def get_weighting_function_satellite(f, layers, surf_emis, mu=1):
    # Compute reflectivity
    refl = 1.0-surf_emis

    down = (get_weighting_function_downwelling(f,
        layers,
        mu)
        .set_index('Zbar').ix[:, ('weighting_function','layer_weight','tx_from_bottom')])
    up = (get_weighting_function_upwelling(f,
        layers,
        mu)
        .set_index('Zbar').ix[:,('weighting_function', 'layer_weight','tx_from_top')])

    # Compose the two weighting functions
    both = down*refl*up.tx_from_top.min() + up

    # Include transmission in data structure for convenience
    both['tx_from_top'] = up.tx_from_top
    both['tx_from_bottom'] = down.tx_from_bottom
    return both


def bt(f, layers, stemp, surf_emis, mu=1.0):
    """
    Compute brightness temperature
    """
    # First get weighting function
    wf = get_weighting_function_satellite(f, layers, surf_emis, mu=mu)
    # Weighted-sum the temperatures and add surface emission
    return wf.tx_from_top.min()*stemp*surf_emis + sum(wf.layer_weight * layers.set_index('Zbar').Tbar)


### All plotting and usage of these functions was documented in an ipython notebook
### This source, the notebook, and report can be found at
### https://github.com/CodaP/aos_weighting_function
