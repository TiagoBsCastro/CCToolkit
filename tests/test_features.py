import numpy as np
import pytest

from cctoolkit.cosmology import CosmologyCalculator
from cctoolkit import baryons
import cctoolkit.utils

# Test CosmologyCalculator initialization and basic attribute computation
def test_cosmology_calculator_initialization():
    params = {
        "H0": 70.0,
        "Ob0": 0.05,
        "Om0": 0.3,
        "sigma8": 0.8,
        "ns": 0.96,
        "TCMB": 2.7255,
        "mnu": 0.06,
        "num_massive_neutrinos": 1,
        "w0": -1.1,
        "wa": 0.8,
    }
    cc = CosmologyCalculator(params)
    # Verify that essential attributes are available
    assert hasattr(cc, "H")
    # Verify that a derived quantity is computed correctly
    omega_m0 = cc.Omega_m(0)
    assert omega_m0 > 0

# Test computation of the halo mass function using both the default and 'castro25' model
def test_halo_mass_function():
    params = {
        "H0": 70.0,
        "Ob0": 0.05,
        "Om0": 0.3,
        "sigma8": 0.8,
        "ns": 0.96,
        "TCMB": 2.7255,
        "mnu": 0.06,
        "num_massive_neutrinos": 1,
        "w0": -1.1,
        "wa": 0.8,
    }
    cc = CosmologyCalculator(params)
    masses = np.logspace(13, 15.5, num=100)
    hmf_default = cc.dndlnM(masses, 0)
    hmf_castro25 = cc.dndlnM(masses, 0, model="castro25")
    # Confirm that the output arrays match the input shape
    assert hmf_default.shape == masses.shape
    assert hmf_castro25.shape == masses.shape
    # Confirm that halo mass function values are positive
    assert np.all(hmf_default > 0)
    assert np.all(hmf_castro25 > 0)

# Test the halo bias functions: PBS prescription and corrected bias
def test_halo_bias():
    params = {
        "H0": 70.0,
        "Ob0": 0.05,
        "Om0": 0.3,
        "sigma8": 0.8,
        "ns": 0.96,
        "TCMB": 2.7255,
        "mnu": 0.06,
        "num_massive_neutrinos": 1,
        "w0": -1.1,
        "wa": 0.8,
    }
    cc = CosmologyCalculator(params)
    masses = np.logspace(13, 15.5, num=100)
    pbs = cc.pbs_bias(masses, 0)
    bias = cc.bias(masses, 0)
    # Ensure the returned arrays match the mass array shape
    assert pbs.shape == masses.shape
    assert bias.shape == masses.shape
    # Verify that the bias values are positive (bias values exceed zero)
    assert np.all(pbs > 0)
    assert np.all(bias > 0)

# Test the conversion of hydro-dynamical masses to dark-matter-only (DMO) virial masses
def test_baryonic_mass_conversion():
    z = 0.0
    params = {"Om0": 0.272,
              "Ob0": 0.272 * 0.168,
              "H0": 70.4,
              "ns": 0.963,
              "mnu": 0,
              "num_massive_neutrinos": 0,
              "sigma8": 0.809}
    cc = CosmologyCalculator(params)
    masses = np.geomspace(1e13, 3e14, num=50)
    # Compute the DMO equivalent mass and corresponding threshold Delta
    M_DMO, Delta = baryons.compute_dmo_mass(masses, z, 0.168)
    # Adjust Delta using the virial overdensity from the utility function
    Delta = Delta * cctoolkit.utils.virial_Delta(cc.Omega_m(z))
    # Convert the threshold to the virial mass using the provided relation
    recovered_masses = [baryons.compute_rec_mass(cc, m, d, z) for m, d in zip(M_DMO, Delta)]
    # Verify that recovered masses have the same number of elements as the input masses
    assert len(recovered_masses) == len(masses)
    # Confirm that all recovered masses are positive
    assert all(mass > 0 for mass in recovered_masses)

# Test the functionality using a tabulated matter power-spectrum
def test_tabulated_power_spectrum():
    masses = np.geomspace(1e13, 1e16, num=50)
    k = np.geomspace(1e-3, 1e1, num=50)
    # Create a dummy power-spectrum that follows a power-law behavior
    Pk = k ** (-2)
    params = {
        "H0": 67.321,
        "Ob0": 0.0494,
        "Om0": 0.3158,
        "sigma8": 0.8102,
        "ns": 0.9661,
        "mnu": 0.0,
        "num_massive_neutrinos": 0
    }
    cc = CosmologyCalculator(params, power_spectrum=[k, Pk])
    dndlnM = cc.dndlnM(masses, 0)
    # Verify that the computed halo mass function has the correct shape and positive values
    assert dndlnM.shape == masses.shape
    assert np.all(dndlnM > 0)

# Test the CosmologyCalculator initialization for a simulation where the radiation contribution is insignificant (low TCMB)
def test_low_TCMB_setting():
    params = {
        "flat": True,
        "H0": 67.321,
        "Om0": 0.3158,
        "Ob0": 0.0494,
        "sigma8": 0.8102,
        "ns": 0.9661,
        "mnu": 0,
        "num_massive_neutrinos": 0,
        "TCMB": 0.5
    }
    cc = CosmologyCalculator(params)
    # Confirm that the calculator initializes properly and computes the growth factor
    assert hasattr(cc, "growth_factor")
