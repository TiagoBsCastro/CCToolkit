"""
hmf.py

This module contains functions and data related to the Halo Mass Function (HMF),
including the multiplicity function and best-fit parameters for different halo finders.
"""

import numpy as np
from scipy.special import gamma

def multiplicity_function(peak_height, dlns_dlnR, Omega_m_z, best_fit_values):
    """
    Calculate the multiplicity function for a given peak height, slope of the power-spectrum,
    and cosmological parameters.

    Parameters:
    -----------
    peak_height : float
        The peak height (nu).
    dlns_dlnR : float
        The logarithmic slope of the variance with respect to the radius.
    Omega_m_z : float
        The matter density parameter at redshift z.
    best_fit_values : dict
        Dictionary containing the best-fit parameters for the chosen halo finder.

    Returns:
    --------
    float
        The value of the multiplicity function.
    
    References
    ----------
    - [arxiv:2208.02174](https://arxiv.org/pdf/2311.01465): "Euclid preparation. XXIV. Calibration of the halo mass function in Λ(ν)CDM cosmologies",
      Castro et al., 2023.
    """
    # Extract best-fit parameters
    a1 = best_fit_values['a1']
    a2 = best_fit_values['a2']
    az = best_fit_values['az']
    p1 = best_fit_values['p1']
    p2 = best_fit_values['p2']
    q1 = best_fit_values['q1']
    q2 = best_fit_values['q2']
    qz = best_fit_values['qz']

    # Compute coefficients aR, qR, a, p, q
    aR = a1 + a2 * (dlns_dlnR + 0.6125) ** 2
    qR = q1 + q2 * (dlns_dlnR + 0.5)
    a = aR * Omega_m_z ** az
    p = p1 + p2 * (dlns_dlnR + 0.5)
    q = qR * Omega_m_z ** qz

    # Amplitude A(p, q)
    A_pq = (2 ** (-1 / 2 - p + q / 2) / np.sqrt(np.pi) * 
            (2 ** p * gamma(q / 2) + gamma(-p + q / 2))) ** (-1)

    # Calculate multiplicity function
    nu = peak_height
    multiplicity = (A_pq * np.sqrt(2 * a * nu ** 2 / np.pi) * 
                    np.exp(-a * nu ** 2 / 2) * 
                    (1 + 1 / (a * nu ** 2) ** p) * 
                    (nu * np.sqrt(a)) ** (q - 1))

    return multiplicity

# Best-fit parameters for different halo finders
best_fit_values_ROCKSTAR = {
    'a1': 0.7962, 'a2': 0.1449, 'az': -0.0658,
    'p1': -0.5612, 'p2': -0.4743,
    'q1': 0.3688, 'q2': -0.2804, 'qz': 0.0251
}

best_fit_values_AHF = {
    'a1': 0.7937, 'a2': 0.1119, 'az': -0.0693,
    'p1': -0.5689, 'p2': -0.4522,
    'q1': 0.3652, 'q2': -0.2628, 'qz': 0.0376
}

best_fit_values_SUBFIND = {
    'a1': 0.7953, 'a2': 0.1667, 'az': -0.0642,
    'p1': -0.6265, 'p2': -0.4907,
    'q1': 0.3215, 'q2': -0.2993, 'qz': 0.0330
}

best_fit_values_VELOCIraptor = {
    'a1': 0.7987, 'a2': 0.1227, 'az': -0.0523,
    'p1': -0.5912, 'p2': -0.4088,
    'q1': 0.3634, 'q2': -0.2732, 'qz': 0.0715
}

