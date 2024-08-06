import cctoolkit
from cctoolkit import baryons
from cctoolkit.cosmology import CosmologyCalculator
import numpy as np
import matplotlib.pyplot as plt

M = np.geomspace(1e13, 1e16, 200)
z = 0.0
calculator_camb = CosmologyCalculator()

M_Delta_dmo, Delta = baryons.compute_dmo_mass(1e14, z, 0.168)
Delta *= cctoolkit.utils.virial_Delta(calculator_camb.Omega_m(z))
baryons.compute_rev_mass(calculator_camb, M_Delta_dmo, Delta, z)
