"""
baryons.py

This module implements the model for the impact of baryonic effects on the masses of clusters
and groups presented in Castro et al. 2024 (https://inspirehep.net/literature/2718844).
"""
import numpy as np
from scipy.optimize import fsolve
from .utils import critical_density0, virial_Delta

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
    m = np.log10(M_vir / 1e10) / 4
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
    m = np.log10(M_vir / 1e10) / 4
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

    delta_f = (0.045 - 0.005 * z) * fb_cosmic  # Baryonic offset parameter
    q = 0.373  # Quasi-adiabatic parameter

    # Step 1: Calculate M∆,dmo using Eq. (5)
    x = (1 - fb_hydro_vir - delta_f) / (1 - fb_cosmic)
    M_Delta_dmo = M_vir_hydro * x

    # Step 3: Calculate the DMO spherical overdensity threshold ∆ using Eq. (4) with respect to the virial Delta
    Delta_over_Delta_vir = x / (1 + q * (1 / x - 1)) ** 3

    return M_Delta_dmo, Delta_over_Delta_vir

def concentration_duffy_et_al_2008 (M, z):

    Avir = 7.85 
    Bvir = -0.081
    Cvir = -0.71
    Mpv  = 2e12

    return Avir * (M/Mpv) ** Bvir * (1+z)**Cvir

def compute_rec_mass(cosmo, M_Delta_dmo, Delta, z, c=concentration_duffy_et_al_2008):

    R_Delta_dmo = np.cbrt(M_Delta_dmo / (4*np.pi*critical_density0*Delta/3))
    Delta_vir = virial_Delta(cosmo.Omega_m(z))
    
    # Auxiliary NFW function
    f = lambda x, _c: np.log(1+_c*x) - _c*x/(1+_c*x)
    # Auxiliary function to be solved
    def fRvir (R): 
        
        Mvir = 4 * np.pi/3 * R**3 * critical_density0 * Delta_vir
        cvir = c(Mvir, z)
        c_Delta = R_Delta_dmo/R * cvir
        
        return M_Delta_dmo * f(R/R_Delta_dmo, c_Delta) / f(1, c_Delta) / (4*np.pi/3*critical_density0*R**3) - Delta_vir
    
    Rvir = fsolve(fRvir, 0.99 * R_Delta_dmo)[0]
    
    return (4*np.pi/3) * Rvir**3 * critical_density0 * Delta_vir