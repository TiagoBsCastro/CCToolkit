"""
baryons.py

This module implements the model for the impact of baryonic effects on the masses of clusters
and groups presented in Castro et al. 2024 (https://inspirehep.net/literature/2718844).
"""
import numpy as np

def baryon_fraction_magneticum(M_vir, z):
    """
    Mean baryon fraction contained inside the virial radius for Magneticum simulations.

    Parameters
    ----------
    M_vir : float
        Virial mass.
    z : float
        Redshift.

    Returns
    -------
    fb_fid : float
        The mean baryon fraction inside the virial radius for Magneticum clusters.
    """
    fb_cosmic = 0.168  # Magneticum baryon fraction
    m = np.log10(M_vir / 1e14)
    gamma = 0.7008 * z - 1.7505
    delta = 0.2
    fb_mean = fb_cosmic * (m ** -gamma) * (1 + m ** (1 / delta)) ** (gamma * delta)
    return fb_mean

def baryon_fraction_tng300(M_vir, z):
    """
    Mean baryon fraction contained inside the virial R for TNG300 simulations.

    Parameters
    ----------
    M_vir : float
        Virial mass.
    z : float
        Redshift.

    Returns
    -------
    fb_fid : float
        The mean baryon fraction inside the virial radius for TNG300 clusters.
    """
    fb_cosmic = 0.1573  # TNG cosmic baryon fraction
    m = np.log10(M_vir / 1e14)
    gamma = 0.4106 * z - 1.4640
    delta = 0.0397 * z ** 2 - 0.04358 * z + 0.03567
    fb_mean = fb_cosmic * (m ** -gamma) * (1 + m ** (1 / delta)) ** (gamma * delta)
    return fb_mean

def compute_dmo_mass(M_vir_hydro, z, fb_cosmic, relation='magneticum'):
    """
    Compute the SO equivalent DMO mass and overdensity Delta (virial units) to
    the virial object definition from hydrodynamic simulations.

    Parameters
    ----------
    M_vir_hydro : float
        Virial mass from hydrodynamic simulations.
    z : float
        Redshift.
    fb_cosmic : float
        Cosmic baryon fraction.
    relation : str, optional
        The baryon fraction relation to use, either 'magneticum' or 'tng300'.
        Default is 'magneticum'.

    Returns
    -------
    M_delta_dmo : float
        The equivalent DMO mass at the Delta overdensity.
    Delta : float
        The overdensity of equivalence in virial units.

    Raises
    ------
    ValueError
        If an invalid relation is specified.
    """
    if relation == 'magneticum':
        fb_hydro_vir = baryon_fraction_magneticum(M_vir_hydro, z)
    elif relation == 'tng300':
        fb_hydro_vir = baryon_fraction_tng300(M_vir_hydro, z)
    else:
        raise ValueError("Invalid relation specified. Use 'magneticum' or 'tng300'.")

    delta_f = 0.045 - 0.005 * z  # Baryonic offset parameter
    q = 0.373  # Quasi-adiabatic parameter

    # Step 1: Calculate M∆,dmo using Eq. (5)
    x = (1 - fb_hydro_vir - delta_f) / (1 - fb_cosmic)
    M_delta_dmo = M_vir_hydro * x

    # Step 3: Calculate the DMO spherical overdensity threshold ∆ using Eq. (4) with respect to the virial Delta
    Delta_over_Delta_vir = x / (1 + q * (1 / x - 1)) ** 3

    return M_delta_dmo, Delta_over_Delta_vir